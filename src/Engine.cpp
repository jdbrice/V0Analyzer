

// STL
#include <iostream>
#include <exception>
#include <stdlib.h> // for atoi

// RooBarb
#include "Logger.h"
#include "LoggerConfig.h"
#include "XmlConfig.h"

#include "TreeAnalyzer.h"
#include "TaskFactory.h"
#include "Engine.h"

using namespace jdb;

// Local
	#include "V0Analyzer.h"

TaskFactory::map_type* TaskFactory::trMap = nullptr;

int main( int argc, char* argv[] ) {

	Logger::setGlobalLogLevel( "all" );
	
	TaskFactory::registerTaskRunner<V0Analyzer>( "V0Analyzer" );

	Engine engine( argc, argv );

	// TaskRunner * tr = TaskFactory::createTaskRunner( "TreeAnalyzer" );
	// tr->run();

	// TaskRunner * tr2 = TaskFactory::createTaskRunner( "TaskRunner" );
	// tr2->run();

	return 0;
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
				// TaskRunner * v0 = new V0Analyzer( config, "V0Analyzer.", fileList, jobPrefix );
				// v0->run();
			} 

		} catch ( exception &e ){
			cout << e.what() << endl;
		}

	}


	return 0;
}
