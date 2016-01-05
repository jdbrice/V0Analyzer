


#include "StHelix.h"
#include "StPhysicalHelixD.h"

void ut_StPhysicalHelix(){

	gSystem->Load( "StHelix_cpp" );

	StPhysicalHelixD test( 0, 0, 0, StThreeVectorD( 0, 0, 0 ), 0 );

}