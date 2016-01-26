#ifndef V0_ANALYZER_H
#define V0_ANALYZER_H 

// JDB 
#include "TreeAnalyzer.h"
#include "XmlConfig.h"
#include "Utils.h"
#include "Logger.h"
#include "CutCollection.h"
using namespace jdb;

// ROOT
#include "TChain.h"

// STL
#include <memory>
#include <vector>
#include <math.h>

// Local
#include "V0PicoDst.h"
#include "Particle.h"


class V0Analyzer : public TreeAnalyzer
{
protected:
	shared_ptr<V0PicoDst> pico;
	
	StThreeVectorD vtx;

	vector<Particle> pos;
	vector<Particle> neg;
	vector<ParticlePair> us;
	vector<ParticlePair*> sig;

	shared_ptr<CutCollection> cc;


	enum class ParticleType { K0S, Lambda, Phi };
	ParticleType plc;

public:
	virtual const char* classname() { return "V0Analyzer"; }
	V0Analyzer() {}
	V0Analyzer( XmlConfig _config, string _nodePath="", int _jobIndex = -1 );
	~V0Analyzer();

	virtual void init( XmlConfig _config, string _nodePath="", int _jobIndex = -1 );

	

	//virtual void make()
protected:

	/* @inherit */
	virtual void preEventLoop();

	/* @inherit */
	virtual void postEventLoop();

	/* @inherit */
	virtual void analyzeEvent( );

	/**
	 * Analyze a track in the current Event
	 * @piTrack	- Track index 
	 */
	virtual void analyzeTrack( Int_t iTrack );
	virtual void makePairs();
	virtual void analyzePairs();

	/* @inherit */
	virtual bool keepEvent();

	/**
	 * Performs track based cuts
	 * @return 	True 	- Passes all selection cuts
	 *          False 	- Fails 1 or more selection cuts
	 */
	virtual bool keepTrack( Int_t iTrack );


	void makeK0SPair( int i, int j);
	void makeLambdaPair( int i, int j);



};


#endif