#ifndef COMMON_H
#define COMMON_H

#include <iostream>

#ifdef PRINT_MODE
    #define PRINT_LOG(x) do { std::cout << "[INFO] " << x << std::endl; } while(0)
#else
    #define PRINT_LOG(x) do { } while(0)
#endif

#ifdef DEBUG_MODE
    #define DEBUG_LOG(x) do { std::cout << "[DEBUG] " << x << std::endl; } while(0)
#else
    #define DEBUG_LOG(x) do { } while(0)
#endif

#endif