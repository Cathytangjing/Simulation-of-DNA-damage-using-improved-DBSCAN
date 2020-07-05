//
// Created by nuccathy on 18-10-10.
//

#ifndef CLUSTERING_CHEMICAL_EXCEPTIONHANDLER_H
#define CLUSTERING_CHEMICAL_EXCEPTIONHANDLER_H
#include "G4VExceptionHandler.hh"

class ExceptionHandler : public G4VExceptionHandler {
 public:
  ExceptionHandler();
  ~ExceptionHandler() override = default;

  G4bool Notify(const char *originOfException,
                const char *exceptionCode,
                G4ExceptionSeverity severity,
                const char *description) override;
 private:

};

#endif //CLUSTERING_CHEMICAL_EXCEPTIONHANDLER_H
