IF( DEFINED ENV{GENFIT_DIR} )
    SET( GENFIT_INCLUDE_DIRS "$ENV{GENFIT_DIR}/include" )
    IF( EXISTS "$ENV{GENFIT_DIR}/lib64" )
        SET( GENFIT_LIBRARIES "$ENV{GENFIT_DIR}/lib64/libgenfit2.so" )
    ELSEIF( EXISTS "$ENV{GENFIT_DIR}/lib" )
        SET( GENFIT_LIBRARIES "$ENV{GENFIT_DIR}/lib/libgenfit2.so" )
    ELSE()
        message( FATAL_ERROR "Genfit library does not exit" )
    ENDIF()
ENDIF()
