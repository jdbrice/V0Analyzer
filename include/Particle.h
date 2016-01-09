#ifndef PARTICLE_H
#define PARTICLE_H


#include "StDcaGeometry.h"
#include "V0PicoDst.h"

#include <TLorentzVector.h>

#include <memory>

class Particle {



public:

	Particle( int id ){
		this->id = id;
		k0Candidate = false;
	}
	~Particle(){}

	int id;
	char charge;

	StThreeVectorD p;
	StPhysicalHelixD helix;
	TLorentzVector lv;

	// filled by pair context
	StThreeVectorD posAtDCA;
	StThreeVectorD pAtDCA;

	bool k0Candidate;


};


class ParticlePair {
public:
	Particle *plc1, *plc2;
	TLorentzVector lv;

	double pointingAngle;
	pair<double,double> pathLengths;
	StThreeVectorD dca;
	StThreeVectorD sVtx;
	StThreeVectorD decLen;
	StThreeVectorD pAtDCA;



};


#endif