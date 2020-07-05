//
// Created by nuccathy on 18-10-10.
//
#include "G4EventManager.hh"
#include "ExceptionHandler.hh"

ExceptionHandler::ExceptionHandler() : G4VExceptionHandler() {

}

G4bool ExceptionHandler::Notify(const char *originOfException,
                                const char *exceptionCode,
                                G4ExceptionSeverity severity,
                                const char *description) {
  G4cout << originOfException << G4endl;

  G4EventManager::GetEventManager()->AbortCurrentEvent();
  return false;
}

