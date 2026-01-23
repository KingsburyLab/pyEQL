SET(Excel_Files
phreeqc.dat
runphreeqc.xls
withcallback.xls
)

SET(Excel_Dir ${EXAMPLES_DIR}/com/excel)

install(FILES ${Excel_Files} DESTINATION ${Excel_Dir})
