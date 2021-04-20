#pragma once
#include <SDKDDKVer.h>
#define WIN_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers

#include <WinSock2.h>
#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#include <glm/glm.hpp>

#include <qapp.h>
#include <sstream>
#include <thread>
#include <mutex>
#include <map>
#include <functional> 
#include <list>
#define IMGUI_DEFINE_MATH_OPERATORS
#include "imgui.h"
#include "imgui_internal.h"
#include "imgui_freetype.h"
#include <filesystem>
#include <cmath>
#include <string>
#include "../timelinefx.h"
