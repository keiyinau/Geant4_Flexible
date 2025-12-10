#include "event.hh"

MyEventAction::MyEventAction(MyRunAction*)
{

}

MyEventAction::~MyEventAction()
{}

void MyEventAction::BeginOfEventAction(const G4Event* aEvent)
{
	G4cout<<">> Begin of Event:"<<aEvent->GetEventID()<<G4endl;
}

void MyEventAction::EndOfEventAction(const G4Event* aEvent)
{

}