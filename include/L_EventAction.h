/*
 * L_EventAction.h
 *
 *  Created on: Oct 2, 2018
 *      Author: vsevolod
 */

#ifndef SRC_L_EVENTACTION_H_
#define SRC_L_EVENTACTION_H_

#include <G4UserEventAction.hh>
#include "globals.hh"
#include "L_RunAction.h"
#include "L_Hit.h"
#include "L_SteppingAction.h"

class G4Event;

class L_RunAction;
class L_SteppingAction;
class L_PrimaryGeneratorAction;

class L_EventAction: public G4UserEventAction {
public:
    L_EventAction(L_RunAction*, L_SteppingAction*);
    virtual ~L_EventAction();
public:
    virtual void  BeginOfEventAction(const G4Event* );
    virtual void    EndOfEventAction(const G4Event* );

    void SetPrimGenerator(L_PrimaryGeneratorAction *gen){_primGenerator = gen;};
private:
    L_RunAction* runAction;
    L_SteppingAction* _steppingAction;
    G4int printModulo;
    G4int theCollectionID;

    L_PrimaryGeneratorAction* _primGenerator;
};

#endif /* SRC_L_EVENTACTION_H_ */
