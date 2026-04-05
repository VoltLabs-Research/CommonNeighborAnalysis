#pragma once

#include <volt/coordination_structures.h>

#include <mutex>

namespace Volt{

inline void ensureCoordinationStructuresInitialized(){
    static std::once_flag initFlag;
    std::call_once(initFlag, []() {
        CoordinationStructures::initializeStructures();
    });
}

}
