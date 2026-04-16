/*
 * beam_sim.cpp - Multi-electrode axisymmetric beam extraction simulation
 *
 * Features: multiple solvers, beam modes (energy/Twiss), solenoid B-field,
 * field diagnostics, convergence tracking, per-electrode currents,
 * energy distribution.
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
#include "epot_gssolver.hpp"
#include "epot_mgsolver.hpp"
#include "particledatabase.hpp"
#include "particlestatistics.hpp"
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
// Multi-electrode support (up to 10)
// =========================================================================
#define MAX_ELECTRODES 10

struct ElectrodeConfig {
    double x_start, x_end, thickness;
    double r_aperture, r_wall, slope, voltage;
};

int             g_num_electrodes = 0;
ElectrodeConfig g_elec[MAX_ELECTRODES];

template<int N>
bool electrode_solid( double x, double y, double z )
{
    if( N >= g_num_electrodes ) return false;
    const ElectrodeConfig &e = g_elec[N];
    if( x < e.x_start || x > e.x_end ) return false;
    if( y < e.r_aperture ) return false;

    // Chamfer: slope > 0 means Pierce-type angled bore opening
    // slope == 0 means rectangular (flat face, straight bore)
    if( e.slope > 1e-6 ) {
        double half = e.thickness * 0.5;
        double x_mid = e.x_start + half;
        if( N == 0 ) {
            // First electrode: downstream chamfer only
            if( x > x_mid ) {
                double r_ch = e.r_aperture + e.slope * ( x - x_mid );
                if( y < r_ch ) return false;
            }
        } else {
            // Other electrodes: symmetric chamfer on both faces
            double r_ch = e.r_aperture + e.slope * fabs( x - x_mid );
            if( y < r_ch ) return false;
        }
    }
    return true;
}

typedef bool (*SolidFunc)( double, double, double );
SolidFunc g_solid_funcs[MAX_ELECTRODES] = {
    electrode_solid<0>, electrode_solid<1>, electrode_solid<2>,
    electrode_solid<3>, electrode_solid<4>, electrode_solid<5>,
    electrode_solid<6>, electrode_solid<7>, electrode_solid<8>,
    electrode_solid<9>
};

// =========================================================================
// Configuration
// =========================================================================
struct SimConfig {
    int    grid_nx      = 241;
    int    n_particles  = 5000;
    double current_density = 600.0;
    double beam_energy  = 5.0;
    double beam_mass    = 1.0;
    double beam_charge  = 1.0;

    // Solver: "bicgstab", "multigrid", "gaussseidel"
    string solver_type  = "bicgstab";
    int    mg_levels    = 4;

    // Beam mode: "energy" or "twiss"
    string beam_mode    = "energy";
    double beam_alpha   = 0.0;
    double beam_beta    = 0.01;
    double beam_emittance = 1e-6;

    // B-field: "none" or "solenoid"
    string bfield_type  = "none";
    double sol_B0       = 0.0;      // Tesla
    double sol_z1       = 0.0;      // start (m)
    double sol_z2       = 0.0;      // end (m)
};

bool read_config( const string &filename, SimConfig &cfg )
{
    ifstream ifs( filename.c_str() );
    if( !ifs.is_open() ) { cerr << "Cannot open: " << filename << endl; return false; }
    string line;
    g_num_electrodes = 0;
    int n_expected = 0;
    bool reading = false;

    while( getline( ifs, line ) ) {
        if( line.empty() || line[0] == '#' ) continue;

        if( line.find("grid_points=")     == 0 ) cfg.grid_nx        = atoi(line.substr(12).c_str());
        else if( line.find("particles=")       == 0 ) cfg.n_particles    = atoi(line.substr(10).c_str());
        else if( line.find("current_density=") == 0 ) cfg.current_density= atof(line.substr(16).c_str());
        else if( line.find("beam_energy=")     == 0 ) cfg.beam_energy    = atof(line.substr(12).c_str());
        else if( line.find("beam_mass=")       == 0 ) cfg.beam_mass      = atof(line.substr(10).c_str());
        else if( line.find("beam_charge=")     == 0 ) cfg.beam_charge    = atof(line.substr(12).c_str());
        else if( line.find("solver=")          == 0 ) cfg.solver_type    = line.substr(7);
        else if( line.find("mg_levels=")       == 0 ) cfg.mg_levels      = atoi(line.substr(10).c_str());
        else if( line.find("beam_mode=")       == 0 ) cfg.beam_mode      = line.substr(10);
        else if( line.find("beam_alpha=")      == 0 ) cfg.beam_alpha     = atof(line.substr(11).c_str());
        else if( line.find("beam_beta=")       == 0 ) cfg.beam_beta      = atof(line.substr(10).c_str());
        else if( line.find("beam_emittance=")  == 0 ) cfg.beam_emittance = atof(line.substr(15).c_str());
        else if( line.find("bfield=")          == 0 ) cfg.bfield_type    = line.substr(7);
        else if( line.find("sol_B0=")          == 0 ) cfg.sol_B0         = atof(line.substr(7).c_str());
        else if( line.find("sol_z1=")          == 0 ) cfg.sol_z1         = atof(line.substr(7).c_str()) * 1e-3;
        else if( line.find("sol_z2=")          == 0 ) cfg.sol_z2         = atof(line.substr(7).c_str()) * 1e-3;
        else if( line.find("electrodes=")      == 0 ) {
            n_expected = atoi( line.substr(11).c_str() );
            reading = true;
        } else if( reading && g_num_electrodes < n_expected
                           && g_num_electrodes < MAX_ELECTRODES ) {
            double pos_mm, r_mm, v_kV, thick_mm, chamfer_deg, wall_mm = 0;
            istringstream iss( line );
            // Format: position aperture voltage thickness chamfer [wall_radius]
            if( !( iss >> pos_mm >> r_mm >> v_kV >> thick_mm >> chamfer_deg ) ) continue;
            iss >> wall_mm;  // optional 6th field
            ElectrodeConfig &e = g_elec[g_num_electrodes];
            e.r_aperture = r_mm * 1e-3;
            e.thickness  = thick_mm * 1e-3;
            if( e.thickness < 0.5e-3 ) e.thickness = 0.5e-3;  // min 0.5 mm
            e.x_start    = pos_mm * 1e-3;
            e.x_end      = e.x_start + e.thickness;
            // Chamfer: 0 = rectangular (flat face), >0 = Pierce-type
            if( chamfer_deg > 85.0 ) chamfer_deg = 85.0;
            e.slope   = ( chamfer_deg > 0.5 ) ? tan( chamfer_deg * M_PI / 180.0 ) : 0.0;
            e.voltage = v_kV * 1e3;
            e.r_wall  = ( wall_mm > r_mm ) ? wall_mm * 1e-3 : 3.0 * e.r_aperture;
            g_num_electrodes++;
        }
    }

    // Sort by position
    for( int i = 0; i < g_num_electrodes-1; i++ )
        for( int j = i+1; j < g_num_electrodes; j++ )
            if( g_elec[j].x_start < g_elec[i].x_start )
                swap( g_elec[i], g_elec[j] );

    // Overlap check
    for( int i = 0; i < g_num_electrodes-1; i++ ) {
        if( g_elec[i].x_end > g_elec[i+1].x_start ) {
            cerr << "ERROR: Electrodes " << i+1 << " and " << i+2 << " overlap!" << endl;
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

    SimConfig cfg;
    if( !read_config( config_file, cfg ) ) { return 1; }

    // ---- Domain ----
    double x_max = 0, r_max = 0, v_min = 0, v_max = 0;
    for( int i = 0; i < g_num_electrodes; i++ ) {
        x_max = max( x_max, g_elec[i].x_end );
        r_max = max( r_max, g_elec[i].r_wall + 1e-3 );
        r_max = max( r_max, g_elec[i].r_aperture * 5.0 );
        v_min = min( v_min, g_elec[i].voltage );
        v_max = max( v_max, g_elec[i].voltage );
    }
    x_max += 2e-3;
    double h = x_max / ( cfg.grid_nx - 1 );
    int grid_ny = (int)round( r_max / h ) + 1;

    cerr << "=== Beam Extraction Simulation ===" << endl;
    cerr << "Solver : " << cfg.solver_type << endl;
    cerr << "Beam   : m=" << cfg.beam_mass << " amu, q=" << cfg.beam_charge << " e" << endl;
    if( cfg.beam_mode == "twiss" )
        cerr << "  Twiss: a=" << cfg.beam_alpha << " b=" << cfg.beam_beta
             << " e=" << cfg.beam_emittance << endl;
    else
        cerr << "  Energy: " << cfg.beam_energy << " eV" << endl;
    if( cfg.bfield_type == "solenoid" )
        cerr << "B-field: solenoid B0=" << cfg.sol_B0 << " T, z="
             << cfg.sol_z1*1e3 << "-" << cfg.sol_z2*1e3 << " mm" << endl;
    cerr << "Grid   : " << cfg.grid_nx << " x " << grid_ny
         << "  h=" << h*1e6 << " um" << endl;
    for( int i = 0; i < g_num_electrodes; i++ )
        cerr << "  Elec#" << i+1 << " pos=" << g_elec[i].x_start*1e3
             << " r=" << g_elec[i].r_aperture*1e3
             << " V=" << g_elec[i].voltage*1e-3 << " kV" << endl;

    // ---- Geometry ----
    Geometry geom( MODE_CYL, Int3D(cfg.grid_nx, grid_ny, 1), Vec3D(0,0,0), h );
    for( int i = 0; i < g_num_electrodes; i++ ) {
        int bid = 7 + i;
        geom.set_solid( bid, new FuncSolid( g_solid_funcs[i] ) );
        geom.set_boundary( bid, Bound(BOUND_DIRICHLET, g_elec[i].voltage) );
    }
    geom.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, g_elec[g_num_electrodes-1].voltage) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN, 0.0) );
    geom.build_mesh();

    // ---- Solver (selectable) ----
    bool is_electron = ( cfg.beam_charge < 0 );
    EpotSolver *solver_ptr = NULL;

    if( cfg.solver_type == "multigrid" ) {
        EpotMGSolver *mg = new EpotMGSolver( geom );
        mg->set_levels( cfg.mg_levels );
        mg->set_mgcycmax( 20 );
        solver_ptr = mg;
        cerr << "Using Multigrid solver (" << cfg.mg_levels << " levels)" << endl;
    } else if( cfg.solver_type == "gaussseidel" ) {
        solver_ptr = new EpotGSSolver( geom );
        cerr << "Using Gauss-Seidel solver" << endl;
    } else {
        solver_ptr = new EpotBiCGSTABSolver( geom );
        cerr << "Using BiCGSTAB solver" << endl;
    }

    InitialPlasma initp( AXIS_X, g_elec[0].thickness * 0.3 );
    if( !is_electron ) {
        solver_ptr->set_initial_plasma( cfg.beam_energy, &initp );
    }

    // ---- Fields ----
    EpotField       epot( geom );
    MeshScalarField scharge( geom );
    EpotEfield      efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    // ---- Magnetic field ----
    MeshVectorField *bfield_ptr = NULL;
    if( cfg.bfield_type == "solenoid" && cfg.sol_B0 != 0 ) {
        bool fout[3] = { true, true, false };  // Bx, Br components, no Btheta
        bfield_ptr = new MeshVectorField( MODE_CYL, fout,
                                          Int3D(cfg.grid_nx, grid_ny, 1),
                                          Vec3D(0,0,0), h );
        for( int i = 0; i < cfg.grid_nx; i++ ) {
            for( int j = 0; j < grid_ny; j++ ) {
                double x = i * h;
                double Bz = ( x >= cfg.sol_z1 && x <= cfg.sol_z2 ) ? cfg.sol_B0 : 0.0;
                bfield_ptr->set( i, j, Vec3D(Bz, 0, 0) );
            }
        }
        cerr << "Solenoid B-field set" << endl;
    } else {
        bfield_ptr = new MeshVectorField;
    }

    // ---- Particles ----
    ParticleDataBaseCyl pdb( geom );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint( true );
    double r_source = g_elec[0].r_wall;

    // ---- Convergence tracking ----
    Convergence conv;
    conv.add_epot( epot );
    conv.add_scharge( scharge );
    Emittance conv_emit;
    conv.add_emittance( 0, conv_emit );

    // ---- Vlasov iteration ----
    int n_iter = 5;
    for( int it = 0; it < n_iter; it++ ) {
        if( it == 1 && !is_electron ) {
            double rhoe = pdb.get_rhosum();
            solver_ptr->set_pexp_plasma( -rhoe, cfg.beam_energy, cfg.beam_energy );
        }

        solver_ptr->solve( epot, scharge );
        efield.recalculate();

        pdb.clear();

        if( cfg.beam_mode == "twiss" ) {
            // Twiss-parameter beam (KV distribution)
            pdb.add_2d_gaussian_beam_with_emittance(
                cfg.n_particles,
                cfg.current_density * r_source,  // total current (A) = J * aperture_width
                cfg.beam_charge, cfg.beam_mass,
                cfg.beam_alpha, cfg.beam_beta, cfg.beam_emittance,
                cfg.beam_energy, 0.0 );
        } else {
            // Energy/temperature beam (default)
            pdb.add_2d_beam_with_energy(
                cfg.n_particles,
                cfg.current_density, cfg.beam_charge, cfg.beam_mass,
                cfg.beam_energy, 0.0, 2.0,
                0.0, 0.0, 0.0, r_source );
        }

        pdb.iterate_trajectories( scharge, efield, *bfield_ptr );

        // Update convergence emittance
        double x_exit = g_elec[g_num_electrodes-1].x_end;
        ParticleDiagPlotter conv_pp( geom, pdb, AXIS_X,
                                     min(x_exit, x_max - 2*h),
                                     PARTICLE_DIAG_PLOT_SCATTER,
                                     DIAG_R, DIAG_RP );
        conv_emit = conv_pp.calculate_emittance();
        conv.evaluate_iteration();

        cerr << "Iteration " << it+1 << "/" << n_iter << " done" << endl;
    }

    // ==================================================================
    // OUTPUT 1 : trajectory plot
    // ==================================================================
    {
        MeshScalarField tdens( geom );
        pdb.build_trajectory_density_field( tdens );
        GeomPlotter gplotter( geom );
        gplotter.set_size( 1200, 600 );
        gplotter.set_font_size( 16 );
        gplotter.set_epot( &epot );
        vector<double> eqlines;
        double vmid = (v_min+v_max)*0.5, vspan = max(v_max-v_min, 1.0);
        for( int k = -4; k <= 4; k++ ) {
            double v = vmid + k * vspan/5.0;
            eqlines.push_back( fabs(v)<1 ? 0.01 : v );
        }
        gplotter.set_eqlines_manual( eqlines );
        gplotter.set_particle_database( &pdb );
        gplotter.set_particle_div( 0 );
        gplotter.set_trajdens( &tdens );
        gplotter.set_fieldgraph_plot( FIELD_TRAJDENS );
        gplotter.fieldgraph()->set_zscale( ZSCALE_RELLOG );
        gplotter.plot_png( outdir + "/trajectory.png" );
    }

    // ==================================================================
    // OUTPUT 1b : domain bounds (for GUI axis labelling)
    // ==================================================================
    {
        ofstream ofs( (outdir + "/domain.csv").c_str() );
        ofs << "x_min_mm,x_max_mm,r_min_mm,r_max_mm\n";
        ofs << 0.0 << "," << x_max*1e3 << "," << 0.0 << "," << r_max*1e3 << "\n";
    }

    // ==================================================================
    // OUTPUT 2 : emittance + phase space + energy
    // ==================================================================
    {
        double x_start = g_elec[0].x_end + h;
        double x_end   = x_max - 3*h;
        int n_samples  = 50;
        double dx = (x_end - x_start) / (n_samples - 1);

        ofstream ofs_e( (outdir + "/emittance_profile.csv").c_str() );
        ofs_e << "x_mm,epsilon,alpha,beta,gamma,r_rms_mm,divergence_mrad,current_A\n";

        ofstream ofs_ps( (outdir + "/phase_space.csv").c_str() );
        ofs_ps << "x_mm,y_mm,yp_mrad\n";

        ofstream ofs_rps( (outdir + "/phase_space_radial.csv").c_str() );
        ofs_rps << "x_mm,r_mm,rp_mrad\n";

        ofstream ofs_ek( (outdir + "/energy_distribution.csv").c_str() );
        ofs_ek << "x_mm,ek_eV\n";

        for( int s = 0; s < n_samples; s++ ) {
            double xp = x_start + s * dx;

            ParticleDiagPlotter pp( geom, pdb, AXIS_X, xp,
                                    PARTICLE_DIAG_PLOT_SCATTER, DIAG_R, DIAG_RP );
            const Emittance &em = pp.calculate_emittance();
            double eps=em.epsilon(), al=em.alpha(), be=em.beta(), ga=em.gamma();
            double r_rms = (eps>0&&be>0) ? sqrt(eps*be) : 0;
            double div   = (eps>0&&ga>0) ? sqrt(eps*ga) : 0;

            ofs_e << fixed << setprecision(6) << xp*1e3 << ","
                  << scientific << setprecision(6) << eps << "," << al << ","
                  << be << "," << ga << ","
                  << fixed << setprecision(6) << r_rms*1e3 << "," << div*1e3 << ","
                  << scientific << setprecision(6) << em.current() << "\n";

            // Y, Y' scatter
            { vector<trajectory_diagnostic_e> d; d.push_back(DIAG_Y); d.push_back(DIAG_YP);
              TrajectoryDiagnosticData td; pdb.trajectories_at_plane(td,AXIS_X,xp,d);
              for( size_t j=0; j<td.traj_size(); j++ )
                  ofs_ps << fixed << setprecision(6) << xp*1e3 << ","
                         << td(j,0)*1e3 << "," << td(j,1)*1e3 << "\n";
            }

            // R, R' scatter (radial phase space, consistent with emittance calc)
            { vector<trajectory_diagnostic_e> d; d.push_back(DIAG_R); d.push_back(DIAG_RP);
              TrajectoryDiagnosticData td; pdb.trajectories_at_plane(td,AXIS_X,xp,d);
              for( size_t j=0; j<td.traj_size(); j++ )
                  ofs_rps << fixed << setprecision(6) << xp*1e3 << ","
                          << td(j,0)*1e3 << "," << td(j,1)*1e3 << "\n";
            }

            // Energy at a few positions
            if( s % 10 == 0 || s == n_samples-1 ) {
                vector<trajectory_diagnostic_e> d; d.push_back(DIAG_EK);
                TrajectoryDiagnosticData td; pdb.trajectories_at_plane(td,AXIS_X,xp,d);
                for( size_t j=0; j<td.traj_size(); j++ )
                    ofs_ek << fixed << setprecision(6) << xp*1e3 << ","
                           << td(j,0) << "\n";
            }
        }
        ofs_e.close(); ofs_ps.close(); ofs_ek.close();
    }

    // ==================================================================
    // OUTPUT 3 : field along axis (potential, E-field, B-field)
    // ==================================================================
    {
        ofstream ofs( (outdir + "/field_along_axis.csv").c_str() );
        ofs << "x_mm,potential_V,Ex_Vm,Bx_T\n";
        int npts = cfg.grid_nx;
        for( int i = 0; i < npts; i++ ) {
            double x = i * h;
            Vec3D pos( x, 0.001, 0 );  // slightly off axis to avoid boundary
            double pot = epot( pos );
            Vec3D ef = efield( pos );
            Vec3D bf(0,0,0);
            if( cfg.bfield_type == "solenoid" && bfield_ptr )
                bf = (*bfield_ptr)( pos );
            ofs << fixed << setprecision(6) << x*1e3 << ","
                << scientific << setprecision(6) << pot << "," << ef[0] << ","
                << bf[0] << "\n";
        }
        ofs.close();
    }

    // ==================================================================
    // OUTPUT 3b : exit-plane particle data (for beam handoff)
    // ==================================================================
    {
        double x_exit = g_elec[g_num_electrodes-1].x_end;
        double xp = min( x_exit + h, x_max - 2*h );
        vector<trajectory_diagnostic_e> diags;
        diags.push_back( DIAG_Y );
        diags.push_back( DIAG_VX );
        diags.push_back( DIAG_VY );
        diags.push_back( DIAG_EK );
        TrajectoryDiagnosticData tdata;
        pdb.trajectories_at_plane( tdata, AXIS_X, xp, diags );

        ofstream ofs( (outdir + "/exit_particles.csv").c_str() );
        ofs << "y_m,vx_ms,vy_ms,vz_ms,ek_eV\n";
        for( size_t j = 0; j < tdata.traj_size(); j++ ) {
            ofs << scientific << setprecision(8)
                << tdata(j,0) << "," << tdata(j,1) << ","
                << tdata(j,2) << "," << 0.0 << ","
                << tdata(j,3) << "\n";
        }
        ofs.close();
        cerr << "Wrote exit_particles.csv (" << tdata.traj_size() << " particles)" << endl;
    }

    // ==================================================================
    // OUTPUT 4 : convergence history
    // ==================================================================
    {
        ofstream ofs( (outdir + "/convergence.csv").c_str() );
        conv.print_history( ofs );
        ofs.close();
    }

    // ==================================================================
    // OUTPUT 5 : per-electrode current
    // ==================================================================
    {
        ofstream ofs( (outdir + "/electrode_currents.csv").c_str() );
        ofs << "electrode,boundary_id,current_A,collisions\n";
        const ParticleStatistics &stats = pdb.get_statistics();
        for( int i = 0; i < g_num_electrodes; i++ ) {
            int bid = 7 + i;
            double cur = 0; uint32_t col = 0;
            if( (uint32_t)bid < stats.number_of_boundaries() ) {
                cur = stats.bound_current( bid );
                col = stats.bound_collisions( bid );
            }
            ofs << i+1 << "," << bid << ","
                << scientific << setprecision(6) << cur << "," << col << "\n";
        }
        // Domain boundaries
        for( int bid = 1; bid <= 4; bid++ ) {
            string name = (bid==1)?"left":(bid==2)?"right":(bid==3)?"axis":"top";
            double cur = 0; uint32_t col = 0;
            if( (uint32_t)bid < stats.number_of_boundaries() ) {
                cur = stats.bound_current( bid );
                col = stats.bound_collisions( bid );
            }
            ofs << name << "," << bid << ","
                << scientific << setprecision(6) << cur << "," << col << "\n";
        }
        ofs.close();
    }

    delete solver_ptr;
    delete bfield_ptr;
    cerr << "Simulation complete." << endl;
    return 0;
}
