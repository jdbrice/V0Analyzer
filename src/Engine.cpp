

// STL
#include <iostream>
#include <exception>
#include <stdlib.h> // for atoi

// RooBarb
#include "Logger.h"
#include "LoggerConfig.h"
#include "XmlConfig.h"
using namespace jdb;

// Local
	#include "V0Analyzer.h"

int main( int argc, char* argv[] ) {

	
	Logger::setGlobalLogLevel( "all" );

	if ( argc >= 2 ){

		try{
			XmlConfig config( argv[ 1 ] );
			//config.report();

			LoggerConfig::setup( &config, "Logger" );

			string fileList = "";
			string jobPrefix = "";

			if ( argc >= 4 ){
				fileList = argv[ 2 ];
				jobPrefix = argv[ 3 ];
			}

			string job = config.getString( "job" );

			if ( "V0Analyzer" == job ){
				V0Analyzer v0( &config, "V0Analyzer.", fileList, jobPrefix );
				v0.make();
			} 

		} catch ( exception &e ){
			cout << e.what() << endl;
		}

	}


	return 0;
}
