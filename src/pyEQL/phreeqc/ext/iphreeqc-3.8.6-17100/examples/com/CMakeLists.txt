add_subdirectory(excel)
add_subdirectory(python)

SET(COM_Files
README.txt
)

SET(COM_Dir ${EXAMPLES_DIR}/com)

install(FILES ${COM_Files} DESTINATION ${COM_Dir})
