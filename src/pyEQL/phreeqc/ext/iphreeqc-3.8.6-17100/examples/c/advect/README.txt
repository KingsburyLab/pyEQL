Example usage

--------------------------------------------------------------------------------------
Windows
--------------------------------------------------------------------------------------
Configure, build, test and install IPhreeqc
  1. cd iphreeqc-3.8.6-17100
  2. mkdir _build
  3. cd _build
  4. cmake -S .. -B . -DCMAKE_INSTALL_PREFIX:PATH=c:/Users/charlton/iphreeqc
  5. cmake --build . --config debug
  6. ctest -C debug
  7. cmake --build . --config release
  8. ctest -C release
  9. cmake --build . --config debug --target install
 10. cmake --build . --config release --target install


Build example:
  1. cd c:\Users\charlton\iphreeqc\examples\c\advect
  2. mkdir _build
  3. cd _build
  4. cmake -DCMAKE_PREFIX_PATH:PATH=c:/Users/charlton/iphreeqc ..
  5. cmake --build . --config release

--------------------------------------------------------------------------------------
Linux/macOS
--------------------------------------------------------------------------------------
Configure, build, test and install IPhreeqc
  1. cd iphreeqc-3.8.6-17100
  2. mkdir _build
  3. cd _build
  4. cmake -S .. -B . -DCMAKE_INSTALL_PREFIX:PATH=/home/charlton/iphreeqc
  5. cmake --build .
  6. ctest
  7. cmake --build . --target install


Build example:
  1. cd /home/charlton/iphreeqc/share/doc/IPhreeqc/examples/c/advect
  2. mkdir _build
  3. cd _build
  4. cmake -DCMAKE_PREFIX_PATH:PATH=/home/charlton/iphreeqc ..
  5. cmake --build .
