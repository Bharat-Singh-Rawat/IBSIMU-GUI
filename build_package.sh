#!/bin/bash
# =======================================================================
# build_package.sh - Build IBSIMU Beam Extraction GUI as a distributable
#
# Usage:
#   ./build_package.sh              # build for current platform
#   ./build_package.sh appimage     # Linux AppImage  (requires appimagetool)
#   ./build_package.sh pyinstaller  # PyInstaller bundle (any platform)
#
# Prerequisites (install once):
#   pip install pyinstaller matplotlib pillow numpy
# =======================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/build"
DIST_DIR="$SCRIPT_DIR/dist"
VENV_DIR="$SCRIPT_DIR/buildenv"

echo "=== IBSIMU GUI Packager ==="
echo "Platform: $(uname -s) $(uname -m)"

# -----------------------------------------------------------------------
# Step 1: Build the C++ simulation binary
# -----------------------------------------------------------------------
echo ""
echo "[1/3] Building C++ simulation backend..."
cd "$SCRIPT_DIR"
make clean 2>/dev/null || true
make
echo "    beam_sim built OK"

# Collect shared libraries that beam_sim needs
collect_libs() {
    local binary="$1"
    local dest="$2"
    mkdir -p "$dest"

    # Copy the binary
    cp "$binary" "$dest/"

    # Copy required shared libs (skip system basics)
    ldd "$binary" 2>/dev/null | while read -r line; do
        lib=$(echo "$line" | awk '{print $3}')
        if [ -f "$lib" ]; then
            case "$lib" in
                /lib/x86_64-linux-gnu/libc.so*) ;;      # skip glibc core
                /lib/x86_64-linux-gnu/libm.so*) ;;
                /lib/x86_64-linux-gnu/libdl.so*) ;;
                /lib/x86_64-linux-gnu/librt.so*) ;;
                /lib/x86_64-linux-gnu/libpthread*) ;;
                /lib64/ld-linux*) ;;
                *)
                    cp -L "$lib" "$dest/" 2>/dev/null || true
                    ;;
            esac
        fi
    done
    echo "    Collected $(ls "$dest"/*.so* 2>/dev/null | wc -l) shared libraries"
}

# -----------------------------------------------------------------------
# Step 2: Bundle with PyInstaller
# -----------------------------------------------------------------------
echo ""
echo "[2/3] Creating PyInstaller bundle..."

# Activate venv if it exists
if [ -d "$VENV_DIR" ]; then
    source "$VENV_DIR/bin/activate"
fi

# Collect beam_sim + its libraries
SIMBIN_DIR="$BUILD_DIR/simbin"
rm -rf "$SIMBIN_DIR"
collect_libs "$SCRIPT_DIR/beam_sim" "$SIMBIN_DIR"

# Create wrapper that sets LD_LIBRARY_PATH for beam_sim
cat > "$BUILD_DIR/beam_sim_wrapper.sh" << 'WRAPPER'
#!/bin/bash
DIR="$(cd "$(dirname "$0")" && pwd)"
export LD_LIBRARY_PATH="$DIR/simbin:$LD_LIBRARY_PATH"
exec "$DIR/simbin/beam_sim" "$@"
WRAPPER
chmod +x "$BUILD_DIR/beam_sim_wrapper.sh"

# Create PyInstaller spec
cat > "$BUILD_DIR/ibsimu_gui.spec" << SPEC
# -*- mode: python ; coding: utf-8 -*-
import os

block_cipher = None
script_dir = '$SCRIPT_DIR'
build_dir = '$BUILD_DIR'

a = Analysis(
    [os.path.join(script_dir, 'beam_gui.py')],
    pathex=[script_dir],
    binaries=[],
    datas=[
        (os.path.join(build_dir, 'simbin'), 'simbin'),
        (os.path.join(build_dir, 'beam_sim_wrapper.sh'), '.'),
    ],
    hiddenimports=[
        'PIL._tkinter_finder',
        'matplotlib.backends.backend_tkagg',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['PyQt5', 'PyQt6', 'PySide2', 'PySide6', 'wx',
              'IPython', 'jupyter', 'notebook', 'scipy', 'pandas'],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='IBSimuGUI',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='IBSimuGUI',
)
SPEC

# Run PyInstaller
cd "$BUILD_DIR"
pyinstaller --distpath "$DIST_DIR" --workpath "$BUILD_DIR/pyinst_work" \
    --clean --noconfirm ibsimu_gui.spec

echo "    PyInstaller bundle created"

# -----------------------------------------------------------------------
# Step 3: Create launch script
# -----------------------------------------------------------------------
echo ""
echo "[3/3] Creating launch script..."

BUNDLE="$DIST_DIR/IBSimuGUI"

# Patch beam_gui.py's SIM_BINARY path for the bundled version
# The bundled app needs to find beam_sim relative to its own location
cat > "$BUNDLE/run_ibsimu.sh" << 'LAUNCH'
#!/bin/bash
# Launch IBSIMU Beam Extraction GUI
DIR="$(cd "$(dirname "$0")" && pwd)"
export LD_LIBRARY_PATH="$DIR/simbin:$LD_LIBRARY_PATH"

# Make beam_sim executable
chmod +x "$DIR/simbin/beam_sim" 2>/dev/null

exec "$DIR/IBSimuGUI" "$@"
LAUNCH
chmod +x "$BUNDLE/run_ibsimu.sh"

# -----------------------------------------------------------------------
# AppImage (optional)
# -----------------------------------------------------------------------
if [ "${1:-}" = "appimage" ]; then
    echo ""
    echo "[Optional] Creating AppImage..."

    APPDIR="$BUILD_DIR/IBSimuGUI.AppDir"
    rm -rf "$APPDIR"
    mkdir -p "$APPDIR/usr/bin" "$APPDIR/usr/lib" "$APPDIR/usr/share/applications" \
             "$APPDIR/usr/share/icons/hicolor/256x256/apps"

    # Copy bundle contents
    cp -r "$BUNDLE"/* "$APPDIR/usr/bin/"

    # Desktop file
    cat > "$APPDIR/usr/share/applications/ibsimu-gui.desktop" << DESKTOP
[Desktop Entry]
Name=IBSIMU Beam Extraction
Exec=IBSimuGUI
Icon=ibsimu
Type=Application
Categories=Science;Physics;
DESKTOP
    cp "$APPDIR/usr/share/applications/ibsimu-gui.desktop" "$APPDIR/"

    # Minimal icon (1x1 pixel PNG as placeholder)
    python3 -c "
from PIL import Image
img = Image.new('RGBA', (256,256), (30,100,200,255))
img.save('$APPDIR/ibsimu.png')
"
    cp "$APPDIR/ibsimu.png" "$APPDIR/usr/share/icons/hicolor/256x256/apps/"

    # AppRun
    cat > "$APPDIR/AppRun" << 'APPRUN'
#!/bin/bash
SELF=$(readlink -f "$0")
HERE=${SELF%/*}
export LD_LIBRARY_PATH="$HERE/usr/bin/simbin:$HERE/usr/lib:$LD_LIBRARY_PATH"
exec "$HERE/usr/bin/IBSimuGUI" "$@"
APPRUN
    chmod +x "$APPDIR/AppRun"

    # Build AppImage (needs appimagetool)
    if command -v appimagetool &>/dev/null; then
        ARCH=x86_64 appimagetool "$APPDIR" "$DIST_DIR/IBSimuGUI-x86_64.AppImage"
        echo "    AppImage created: $DIST_DIR/IBSimuGUI-x86_64.AppImage"
    else
        echo "    appimagetool not found. Install from:"
        echo "    https://github.com/AppImage/AppImageKit/releases"
        echo "    AppDir prepared at: $APPDIR"
    fi
fi

# -----------------------------------------------------------------------
echo ""
echo "=== Build Complete ==="
echo ""
echo "Distribution at: $DIST_DIR/IBSimuGUI/"
echo ""
echo "To run:  $DIST_DIR/IBSimuGUI/run_ibsimu.sh"
echo ""
echo "To distribute: zip/tar the IBSimuGUI/ folder."
echo "  cd $DIST_DIR && tar czf IBSimuGUI-linux-x86_64.tar.gz IBSimuGUI/"
echo ""
