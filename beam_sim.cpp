/*
 * beam_sim.cpp - Multi-electrode axisymmetric beam extraction simulation
 *
 * Reads electrode definitions from a config file, runs Vlasov iteration,
 * and outputs trajectory plot + emittance/phase-space data for the GUI.
 */

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "epot_bicgstabsolver.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "convergence.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"
#include "geomplotter.hpp"
#include "trajectorydiagnostics.hpp"
#include "config.h"

using namespace std;

// =========================================================================
// Multi-electrode support (up to 10 electrodes)
// =========================================================================
#define MAX_ELECTRODES 10

struct ElectrodeConfig {
    double x_start;      // front face position (m)
    double x_end;        // back face position (m)
    double thickness;    // plate thickness (m)
    double r_aperture;   // bore aperture radius (m)
    double r_wall;       // back-wall radius for first electrode (m)
    double slope;        // chamfer slope = tan(chamfer_angle)
    double voltage;      // voltage (V), sign included
};

int              g_num_electrodes = 0;
ElectrodeConfig  g_elec[MAX_ELECTRODES];

// ---------------------------------------------------------------------------
// Generic electrode solid (template index selects which electrode)
// ---------------------------------------------------------------------------
template<int N>
bool electrode_solid( double x, double y, double z )
{
    if( N >= g_num_electrodes ) return false;
    const ElectrodeConfig &e = g_elec[N];

    if( x < e.x_start || x > e.x_end ) return false;
    if( y < e.r_aperture )              return false;

    double half  = e.thickness * 0.5;
    double x_mid = e.x_start + half;

    if( N == 0 ) {
        // First electrode: back-wall on upstream face, chamfer on downstream
        double x_wall_end = e.x_start + e.thickness * 0.25;
        if( x < x_wall_end && y < e.r_wall ) return false;

        if( x > x_mid ) {
            double r_ch = e.r_aperture + e.slope * ( x - x_mid );
            if( y < r_ch ) return false;
        }
    } else {
        // Other electrodes: symmetric chamfer on both faces
        double dist_from_mid = fabs( x - x_mid );
        double r_ch = e.r_aperture + e.slope * dist_from_mid;
        if( y < r_ch ) return false;
    }
    return true;
}

// Function-pointer table (template instantiations)
typedef bool (*SolidFunc)( double, double, double );
SolidFunc g_solid_funcs[MAX_ELECTRODES] = {
    electrode_solid<0>, electrode_solid<1>, electrode_solid<2>,
    electrode_solid<3>, electrode_solid<4>, electrode_solid<5>,
    electrode_solid<6>, electrode_solid<7>, electrode_solid<8>,
    electrode_solid<9>
};

// =========================================================================
// Config-file reader
// =========================================================================
// Format:
//   grid_points=241
//   particles=5000
//   current_density=600.0
//   beam_energy=5.0
//   electrodes=2
//   0.0  0.5  0.0   63.4      (position_mm  aperture_mm  voltage_kV  chamfer_deg)
//   10.0 1.5  -8.0  45.0
// =========================================================================
bool read_config( const string &filename,
                  int &grid_nx, int &n_particles,
                  double &current_density, double &beam_energy,
                  double &beam_mass, double &beam_charge )
{
    ifstream ifs( filename.c_str() );
    if( !ifs.is_open() ) {
        cerr << "Cannot open config: " << filename << endl;
        return false;
    }
    string line;
    g_num_electrodes = 0;
    int n_expected = 0;
    bool reading = false;

    while( getline( ifs, line ) ) {
        if( line.empty() || line[0] == '#' ) continue;

        if( line.find("grid_points=")    == 0 ) grid_nx        = atoi( line.substr(12).c_str() );
        else if( line.find("particles=")      == 0 ) n_particles    = atoi( line.substr(10).c_str() );
        else if( line.find("current_density=")== 0 ) current_density= atof( line.substr(16).c_str() );
        else if( line.find("beam_energy=")    == 0 ) beam_energy    = atof( line.substr(12).c_str() );
        else if( line.find("beam_mass=")      == 0 ) beam_mass      = atof( line.substr(10).c_str() );
        else if( line.find("beam_charge=")    == 0 ) beam_charge    = atof( line.substr(12).c_str() );
        else if( line.find("electrodes=")     == 0 ) {
            n_expected = atoi( line.substr(11).c_str() );
            reading = true;
        } else if( reading && g_num_electrodes < n_expected
                           && g_num_electrodes < MAX_ELECTRODES ) {
            double pos_mm, r_mm, v_kV, chamfer_deg;
            istringstream iss( line );
            if( !( iss >> pos_mm >> r_mm >> v_kV >> chamfer_deg ) ) {
                cerr << "Bad electrode line: " << line << endl;
                continue;
            }
            ElectrodeConfig &e = g_elec[g_num_electrodes];
            e.r_aperture = r_mm  * 1e-3;
            e.thickness  = max( 2.0e-3, 4.0 * e.r_aperture );
            e.x_start    = pos_mm * 1e-3;
            e.x_end      = e.x_start + e.thickness;
            // Clamp chamfer to avoid zero/negative slopes
            if( chamfer_deg < 5.0 )  chamfer_deg = 5.0;
            if( chamfer_deg > 85.0 ) chamfer_deg = 85.0;
            e.slope      = tan( chamfer_deg * M_PI / 180.0 );
            e.voltage    = v_kV * 1e3;       // sign from user (e.g. -8 kV → -8000 V)
            e.r_wall     = 3.0 * e.r_aperture;
            g_num_electrodes++;
        }
    }

    // Sort electrodes by position
    for( int i = 0; i < g_num_electrodes - 1; i++ )
        for( int j = i + 1; j < g_num_electrodes; j++ )
            if( g_elec[j].x_start < g_elec[i].x_start )
                swap( g_elec[i], g_elec[j] );

    // Check for overlapping electrodes
    for( int i = 0; i < g_num_electrodes - 1; i++ ) {
        if( g_elec[i].x_end > g_elec[i+1].x_start ) {
            cerr << "ERROR: Electrodes " << i+1 << " and " << i+2
                 << " overlap! Electrode " << i+1 << " ends at "
                 << g_elec[i].x_end*1e3 << " mm but electrode " << i+2
                 << " starts at " << g_elec[i+1].x_start*1e3 << " mm."
                 << endl;
            cerr << "Increase the distance between electrodes or reduce"
                    " aperture sizes." << endl;
            return false;
        }
    }

    return g_num_electrodes > 0;
}

// =========================================================================
// Main
// =========================================================================
int main( int argc, char **argv )
{
    string config_file = "beam_config.txt";
    string outdir      = ".";

    for( int i = 1; i < argc; i++ ) {
        if( !strcmp(argv[i],"--config") && i+1<argc ) config_file = argv[++i];
        else if( !strcmp(argv[i],"--outdir") && i+1<argc ) outdir = argv[++i];
    }

    // Defaults (overridden by config)
    int    grid_nx         = 241;
    int    n_particles     = 5000;
    double current_density = 600.0;
    double beam_energy     = 5.0;
    double beam_mass       = 1.0;    // amu (1 = proton)
    double beam_charge     = 1.0;    // multiples of e (1 = singly charged)

    if( !read_config( config_file, grid_nx, n_particles,
                      current_density, beam_energy,
                      beam_mass, beam_charge ) ) {
        cerr << "Failed to read config file: " << config_file << endl;
        return 1;
    }

    // ---- Compute domain size ----
    double x_max = 0, r_max = 0;
    double v_min = 0, v_max = 0;
    for( int i = 0; i < g_num_electrodes; i++ ) {
        x_max = max( x_max, g_elec[i].x_end );
        r_max = max( r_max, g_elec[i].r_wall  + 1.0e-3 );
        r_max = max( r_max, g_elec[i].r_aperture * 5.0 );
        v_min = min( v_min, g_elec[i].voltage );
        v_max = max( v_max, g_elec[i].voltage );
    }
    x_max += 2.0e-3;   // padding

    double h      = x_max / ( grid_nx - 1 );
    int    grid_ny = (int)round( r_max / h ) + 1;

    cerr << "=== Multi-Electrode Beam Extraction ===" << endl;
    cerr << "Electrodes: " << g_num_electrodes << endl;
    for( int i = 0; i < g_num_electrodes; i++ )
        cerr << "  #" << i+1
             << "  pos=" << g_elec[i].x_start*1e3 << " mm"
             << "  r=" << g_elec[i].r_aperture*1e3 << " mm"
             << "  V=" << g_elec[i].voltage*1e-3 << " kV"
             << "  chamfer=" << atan(g_elec[i].slope)*180.0/M_PI << " deg"
             << endl;
    cerr << "Beam   : m=" << beam_mass << " amu, q=" << beam_charge
         << " e, E=" << beam_energy << " eV" << endl;
    cerr << "Grid   : " << grid_nx << " x " << grid_ny
         << "  h=" << h*1e6 << " um" << endl;
    cerr << "Domain : " << x_max*1e3 << " x " << r_max*1e3 << " mm" << endl;

    // ---- Build geometry (cylindrical symmetry) ----
    Geometry geom( MODE_CYL, Int3D(grid_nx, grid_ny, 1),
                   Vec3D(0,0,0), h );

    for( int i = 0; i < g_num_electrodes; i++ ) {
        int bid = 7 + i;
        Solid *s = new FuncSolid( g_solid_funcs[i] );
        geom.set_solid( bid, s );
        geom.set_boundary( bid, Bound(BOUND_DIRICHLET, g_elec[i].voltage) );
    }

    geom.set_boundary( 1, Bound(BOUND_NEUMANN,   0.0) );   // left
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,
                                g_elec[g_num_electrodes-1].voltage) ); // right
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,   0.0) );   // axis
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,   0.0) );   // top
    geom.build_mesh();

    // ---- Solver ----
    EpotBiCGSTABSolver solver( geom );
    bool is_electron = ( beam_charge < 0 );

    // Plasma expansion model only valid for positive ions
    if( !is_electron ) {
        InitialPlasma initp( AXIS_X, g_elec[0].thickness * 0.3 );
        solver.set_initial_plasma( beam_energy, &initp );
    }

    EpotField       epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield      efield( epot );

    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBaseCyl pdb( geom );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );

    double r_source = g_elec[0].r_wall;

    // ---- Vlasov iteration ----
    int n_iter = 5;
    for( int it = 0; it < n_iter; it++ ) {

        // Plasma expansion model for ions only (not electrons)
        if( it == 1 && !is_electron ) {
            double rhoe = pdb.get_rhosum();
            solver.set_pexp_plasma( -rhoe, beam_energy, beam_energy );
        }

        solver.solve( epot, scharge );
        efield.recalculate();

        pdb.clear();
        pdb.add_2d_beam_with_energy( n_particles,
                                     current_density, beam_charge, beam_mass,
                                     beam_energy, 0.0, 2.0,
                                     0.0, 0.0,
                                     0.0, r_source );
        pdb.iterate_trajectories( scharge, efield, bfield );
        cerr << "Iteration " << it+1 << "/" << n_iter << " done" << endl;
    }

    // ==================================================================
    // OUTPUT 1 : trajectory geometry plot
    // ==================================================================
    {
        MeshScalarField tdens( geom );
        pdb.build_trajectory_density_field( tdens );

        GeomPlotter gplotter( geom );
        gplotter.set_size( 1200, 600 );
        gplotter.set_font_size( 16 );
        gplotter.set_epot( &epot );

        vector<double> eqlines;
        double vmid = ( v_min + v_max ) * 0.5;
        double vspan = v_max - v_min;
        if( vspan < 1.0 ) vspan = 1.0;
        for( int k = -4; k <= 4; k++ ) {
            double v = vmid + k * vspan / 5.0;
            eqlines.push_back( fabs(v) < 1.0 ? 0.01 : v );
        }
        gplotter.set_eqlines_manual( eqlines );
        gplotter.set_particle_database( &pdb );
        gplotter.set_particle_div( 0 );
        gplotter.set_trajdens( &tdens );
        gplotter.set_fieldgraph_plot( FIELD_TRAJDENS );
        gplotter.fieldgraph()->set_zscale( ZSCALE_RELLOG );

        string fname = outdir + "/trajectory.png";
        gplotter.plot_png( fname );
        cerr << "Wrote " << fname << endl;
    }

    // ==================================================================
    // OUTPUT 2 : emittance profile + phase-space scatter (Y, Y')
    // ==================================================================
    {
        double x_start = g_elec[0].x_end + h;
        double x_end   = x_max - 3.0 * h;
        int    n_samples = 50;
        double dx = ( x_end - x_start ) / ( n_samples - 1 );

        // -- emittance profile --
        string efile = outdir + "/emittance_profile.csv";
        ofstream ofs_emit( efile.c_str() );
        ofs_emit << "x_mm,epsilon,alpha,beta,gamma,"
                    "r_rms_mm,divergence_mrad,current_A\n";

        // -- phase-space scatter (Y, Y') --
        string pfile = outdir + "/phase_space.csv";
        ofstream ofs_ps( pfile.c_str() );
        ofs_ps << "x_mm,y_mm,yp_mrad\n";

        for( int s = 0; s < n_samples; s++ ) {
            double xp = x_start + s * dx;

            // Emittance via R, R' (proper for cylindrical)
            ParticleDiagPlotter pp( geom, pdb, AXIS_X, xp,
                                    PARTICLE_DIAG_PLOT_SCATTER,
                                    DIAG_R, DIAG_RP );
            const Emittance &em = pp.calculate_emittance();

            double eps     = em.epsilon();
            double alpha_t = em.alpha();
            double beta_t  = em.beta();
            double gamma_t = em.gamma();
            double r_rms   = ( eps > 0 && beta_t > 0 ) ? sqrt(eps*beta_t) : 0;
            double div_rms = ( eps > 0 && gamma_t > 0 ) ? sqrt(eps*gamma_t) : 0;
            double curr    = em.current();

            ofs_emit << fixed      << setprecision(6) << xp*1e3 << ","
                     << scientific << setprecision(6)
                     << eps << "," << alpha_t << "," << beta_t << ","
                     << gamma_t << ","
                     << fixed << setprecision(6)
                     << r_rms*1e3 << "," << div_rms*1e3 << ","
                     << scientific << setprecision(6) << curr << "\n";

            // Per-particle scatter using Y, Y' (gives full ± distribution)
            vector<trajectory_diagnostic_e> diags;
            diags.push_back( DIAG_Y );
            diags.push_back( DIAG_YP );
            TrajectoryDiagnosticData tdata;
            pdb.trajectories_at_plane( tdata, AXIS_X, xp, diags );

            for( size_t j = 0; j < tdata.traj_size(); j++ ) {
                ofs_ps << fixed << setprecision(6)
                       << xp*1e3 << ","
                       << tdata(j,0)*1e3 << ","
                       << tdata(j,1)*1e3 << "\n";
            }
        }
        ofs_emit.close();
        ofs_ps.close();
        cerr << "Wrote " << efile << endl;
        cerr << "Wrote " << pfile << endl;
    }

    cerr << "Simulation complete." << endl;
    return 0;
}
