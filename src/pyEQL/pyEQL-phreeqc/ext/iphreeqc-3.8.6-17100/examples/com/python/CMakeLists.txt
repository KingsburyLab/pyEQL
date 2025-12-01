SET(Python_Files
Gypsum.py
parallel_advect.py
phreeqc.dat
pitzer.dat
wateq4f.dat
)

SET(Python_Dir ${EXAMPLES_DIR}/com/python)

install(FILES ${Python_Files} DESTINATION ${Python_Dir})
