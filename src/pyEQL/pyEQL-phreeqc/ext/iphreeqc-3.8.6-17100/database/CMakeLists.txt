set(phreeqc_DATABASE
  Amm.dat
  ColdChem.dat
  Concrete_PHR.dat
  Concrete_PZ.dat
  core10.dat
  frezchem.dat
  iso.dat
  Kinec_v3.dat
  Kinec.v2.dat
  llnl.dat
  minteq.dat
  minteq.v4.dat
  phreeqc_rates.dat
  PHREEQC_ThermoddemV1.10_15Dec2020.dat
  phreeqc.dat
  pitzer.dat
  sit.dat
  Tipping_Hurley.dat
  wateq4f.dat
)

# for mytest tests
foreach(db ${phreeqc_DATABASE})
  configure_file(${db} ${db} COPYONLY)
endforeach()

if(WIN32)
  install(FILES ${phreeqc_DATABASE} DESTINATION database)
else()
  install(FILES ${phreeqc_DATABASE} DESTINATION ${CMAKE_INSTALL_DOCDIR}/database)
endif()
