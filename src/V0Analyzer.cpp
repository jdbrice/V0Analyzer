#include "V0Analyzer.h"
#include "PhysicalConstants.h"

const float MASS_PI = 0.13957;
const float MASS_P = 0.93827231;

V0Analyzer::V0Analyzer( XmlConfig config, string np, int _jobIndex ) : TreeAnalyzer( config, np, _jobIndex  ){
	DEBUG( classname(), "" );

	init( config, np, _jobIndex );
}



V0Analyzer::~V0Analyzer(){
	DEBUG( classname(), "" );
}

void V0Analyzer::init( XmlConfig config, string nodePath, int _jobIndex ) {
	DEBUG( classname(), "" );
	TreeAnalyzer::init( config, nodePath, _jobIndex );

	Logger::setGlobalLogLevel( config[ nodePath + ".Logger:globalLogLevel" ] );
	Logger::setGlobalColor( true );

	pico = shared_ptr<V0PicoDst>( new V0PicoDst( chain ) );
	cc = shared_ptr<CutCollection>( new CutCollection( config, nodePath + "CutCollection" ) );

	plc = ParticleType::K0S;

	if ( "Lambda" == config[ nodePath + ":plc" ] || "L" == config[ nodePath + ":plc" ] || "l" == config[ nodePath + ":plc" ] )
		plc = ParticleType::Lambda;
}



void V0Analyzer::preEventLoop(){
	DEBUG( classname(), "" );
	TreeAnalyzer::preEventLoop();

	book->cd();
	
}

void V0Analyzer::postEventLoop(){
	DEBUG( classname(), "" );
}


void V0Analyzer::analyzeEvent(){
	DEBUG( classname(), "" );

	pos.clear();
	neg.clear();

	
	vtx = StThreeVectorD(pico->vtx_mX1, pico->vtx_mX2, pico->vtx_mX3 );

	
	// Number of tracks
	Int_t nTracks = pico->nTracks;

	for ( Int_t iTrack = 0; iTrack < nTracks; iTrack ++ ){

		if ( !keepTrack( iTrack ) )
			continue;
		
		analyzeTrack( iTrack );			
	}

	makePairs();
	analyzePairs();

	for ( Particle &p : pos ){
		if ( !p.k0Candidate ) continue;
		book->fill( "piBeta", p.p.mag(), 1.0/pico->tracks_mBeta[ p.id ] );
		book->fill( "piDedx", p.p.mag(), pico->tracks_mDedx[ p.id ] );
	}
	for ( Particle &p : neg ){
		if ( !p.k0Candidate ) continue;
		book->fill( "piBeta", p.p.mag(), 1.0/pico->tracks_mBeta[ p.id ] );
		book->fill( "piDedx", p.p.mag(), pico->tracks_mDedx[ p.id ] );			
	}
}


void V0Analyzer::analyzePairs(){
	TRACE( classname(), "" );
	CutCollection cuts = *cc;

	ParticlePair* pp = nullptr;

	int n = us.size();
	TRACE( classname(), n << " == us.size()" );
	for ( int iPair = 0; iPair < n; iPair++ ){

		pp = &us[ iPair ];
		Particle * p1 = pp->plc1;
		Particle * p2 = pp->plc2;

		pp->plc1->lv.SetPtEtaPhiM( pp->plc1->p.perp(), pp->plc1->p.pseudoRapidity(), pp->plc1->p.phi(), pp->plc1->mass );
		pp->plc2->lv.SetPtEtaPhiM( pp->plc2->p.perp(), pp->plc2->p.pseudoRapidity(), pp->plc2->p.phi(), pp->plc2->mass );

		pp->lv = pp->plc1->lv + pp->plc2->lv;

		book->fill( "invMassUS", pp->lv.M() );

		if ( !cuts[ "invMass" ]->inInclusiveRange( pp->lv.M() ) ) continue;

		pp->pathLengths = p1->helix.pathLengths( pp->plc2->helix );
		p1->posAtDCA = p1->helix.at( pp->pathLengths.first );
		p2->posAtDCA = p2->helix.at( pp->pathLengths.second );


		pp->dca = p1->posAtDCA - p2->posAtDCA;
		pp->sVtx = p2->posAtDCA + pp->dca * 0.5; 	// secondary vertex position
		pp->decLen = pp->sVtx - vtx;				// decay length vector

		

		
		book->fill( "dca", pp->lv.M(), pp->dca.mag() );
		book->fill( "decLen", pp->lv.M(), pp->decLen.mag() );



		if ( pp->decLen.mag() < cuts[ "decLen" ]->min ) continue;
		if ( pp->dca.mag() > cuts[ "dcaMag" ]->max ) continue;
		



		p1->pAtDCA = p1->helix.momentumAt( pp->pathLengths.first, pico->bField / 1000.0 * kilogauss );
		p2->pAtDCA = p2->helix.momentumAt( pp->pathLengths.second, pico->bField / 1000.0 * kilogauss );

		pp->pAtDCA = p1->pAtDCA + p2->pAtDCA;

		pp->pointingAngle = fabs( pp->pAtDCA.angle( pp->decLen ) );

		book->fill( "sig_point", pp->lv.M(), pp->pointingAngle );

		if ( pp->pointingAngle > cuts[ "pointingAngle" ]->max ) continue;

		

		// reset the p1 and p2 lv to use the mom at DCA
		p1->lv.SetPtEtaPhiM( p1->pAtDCA.perp(), p1->pAtDCA.pseudoRapidity(), p1->pAtDCA.phi(), pp->plc1->mass );
		p2->lv.SetPtEtaPhiM( p2->pAtDCA.perp(), p2->pAtDCA.pseudoRapidity(), p2->pAtDCA.phi(), pp->plc2->mass );
		pp->lv = p1->lv + p2->lv;

		book->fill( "sig", pp->lv.M() );

		if ( cuts[ "candInvMass" ]->inInclusiveRange( pp->lv.M() ) ){
			p1->k0Candidate = true;
			p2->k0Candidate = true;
		}
	}
}
void V0Analyzer::analyzeTrack( int iTrack ){
}

bool V0Analyzer::keepEvent( ){

	// none for now;
	return true;
}


void V0Analyzer::makePairs(){
	TRACE( classname(), "" );

	CutCollection cuts = *cc;

	us.clear();

	int npos = pos.size();
	int nneg = neg.size();

	TRACE( classname(), "#Pos = " << npos );
	TRACE( classname(), "#Neg = " << nneg );

	for ( int i = 0; i < npos; i++ ){
		for( int j = 0; j < nneg; j++ ){

			if ( ParticleType::K0S == plc ){
				makeK0SPair( i, j );
			} else if ( ParticleType::Lambda == plc ){
				makeLambdaPair( i, j );
			}			

		} // j
	} // i
}


void V0Analyzer::makeK0SPair( int i, int j ){

	ParticlePair pp;

	pp.plc1 = &pos[i];
	pp.plc2 = &neg[j];

	pp.plc1->mass = MASS_PI;
	pp.plc2->mass = MASS_PI;

	us.push_back( pp );  
}

void V0Analyzer::makeLambdaPair( int i, int j){

	ParticlePair pp1;

	pp1.plc1 = &pos[i];
	pp1.plc2 = &neg[j];

	pp1.plc1->mass = MASS_P;
	pp1.plc2->mass = MASS_PI;

	us.push_back( pp1 );

	ParticlePair pp2;

	pp2.plc1 = &pos[i];
	pp2.plc2 = &neg[j];

	pp2.plc1->mass = MASS_PI;
	pp2.plc2->mass = MASS_P;

	us.push_back( pp2 );

}


bool V0Analyzer::keepTrack( int iTrack ){
	

	DEBUG( classname(), "" );

	CutCollection cuts = *cc;
	Particle plc( iTrack );
	float fitRatio = (float)pico->tracks_mNHitsFit[ iTrack ] / pico->tracks_mNHitsPoss[ iTrack ];
	// book->fill( "pre_nHitsFit", pico->tracks_mNHitsFit[ iTrack ] );
	// book->fill( "pre_nHitsDedx", pico->tracks_mNHitsDedx[ iTrack ] );
	// book->fill( "pre_nHitsPoss", pico->tracks_mNHitsPoss[ iTrack ] );
	// book->fill( "pre_nHitsRatio", fitRatio );	

	plc.charge = pico->tracks_mNHitsFit[iTrack] > 0 ? 1 : -1;

	if ( abs(pico->tracks_mNHitsFit[ iTrack ]) < 16 )
		return false;
	if ( abs(fitRatio) < 0.52 )
		return false;

	float px = pico->tracks_mMomentum_fX[iTrack];
	float py = pico->tracks_mMomentum_fY[iTrack];
	float pz = pico->tracks_mMomentum_fZ[iTrack];

	if ( sqrt( px*px + py*py + pz*pz ) < 0  )
		return false;
	if ( sqrt( px*px + py*py ) < 0.1 )
		return false;

	book->fill( "charge", plc.charge );


	StDcaGeometry dcaGeom;
	dcaGeom.set(pico->tracks_mDcaParams[iTrack], pico->tracks_mDcaParams[iTrack]);

	// global momentum
	StThreeVectorD  otherP = StThreeVectorD( pico->tracks_mMomentum_fX[iTrack], pico->tracks_mMomentum_fY[iTrack], pico->tracks_mMomentum_fZ[iTrack] );

	plc.helix = dcaGeom.helix();
    plc.p = plc.helix.momentum( pico->bField / 1000.0 * kilogauss  );

    // global momentum versus helix momentum - should they be the same?
    book->fill( "pvp", plc.p.magnitude(), otherP.magnitude() );


    if ( plc.p.magnitude() <= 0 ) return false;
    if ( plc.p.perp() < 0.1 ) return false;
    if( fabs( plc.p.pseudoRapidity() ) > 1.0 ) return false;

    
    // DCA cut
    if( plc.helix.distance( vtx ) < cuts[ "dcaDaughterMag" ]->min / pow( plc.p.mag(), 1.5 )   ) return false;

    book->fill( "origin", dcaGeom.origin().magnitude() );
    book->fill( "vd_p", plc.p.mag(), plc.helix.distance( vtx ) );

	// book->fill( "nHitsFit", pico->tracks_mNHitsFit[ iTrack ] );
	// book->fill( "nHitsDedx", pico->tracks_mNHitsDedx[ iTrack ] );
	// book->fill( "nHitsPoss", pico->tracks_mNHitsPoss[ iTrack ] );
	// book->fill( "nHitsRatio", fitRatio );	

	
	if ( plc.charge > 0 )
		pos.push_back( plc );
	else if ( plc.charge < 0 )
		neg.push_back( plc );

	return true;
}
