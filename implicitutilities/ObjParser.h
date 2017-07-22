#ifndef OBJ_PARSER_H
#define OBJ_PARSER_H

#include "MathDefines.h"

namespace ObjParser
{
  void loadMesh( const std::string& file_name, Matrix3Xsc& vertices, Matrix2Xuc& edges, Matrix3Xuc& faces );
}

#endif
