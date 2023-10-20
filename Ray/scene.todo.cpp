#include <cmath>
#include <Util/exceptions.h>
#include "scene.h"

using namespace Ray;
using namespace Util;
using namespace std;

///////////
// Scene //
///////////
Point3D Scene::Reflect( Point3D v , Point3D n )
{
	//////////////////
	// Reflect here //
	//////////////////
	WARN_ONCE( "method undefined" );
	return Point3D();
}

bool Scene::Refract( Point3D v , Point3D n , double ir , Point3D& refract )
{
	//////////////////
	// Refract here //
	//////////////////
	WARN_ONCE( "method undefined" );
	return false;
}

// called directly by the scene, which calls processFirstIntersection's
Point3D Scene::getColor( Ray3D ray , int rDepth , Point3D cLimit , unsigned int lightSamples , unsigned int tIdx )
{
	Point3D color;
	RayTracingStats::IncrementRayNum();
	ShapeProcessingInfo spInfo;
	RayIntersectionFilter rFilter = []( double ){ return true; };
	// sets color of intersected ray
	RayIntersectionKernel rKernel = [&]( const ShapeProcessingInfo &spInfo , const RayShapeIntersectionInfo &_iInfo ) {
		/////////////////////////////////////////////////////////
		// Create the computational kernel that gets the color //
		/////////////////////////////////////////////////////////

		color += spInfo.material->emissive;

		// loop over light sources in scene
		std::vector< Light * > lights = _globalData.lights;
		for (int i = 0; i < lights.size(); i++) {
			Light *light = lights[i];
			color += light->getAmbient(ray, _iInfo, *spInfo.material);
			if (light->isInShadow(_iInfo, *this, tIdx)) continue;
			color += light->getDiffuse(ray, _iInfo, *spInfo.material);
			color += light->getSpecular(ray, _iInfo, *spInfo.material);
		}

		// Point3D emissive = spInfo.material->emissive;
		// Point3D ambient = spInfo.material->ambient;
		// Point3D diffuse = spInfo.material->diffuse;
		// Point3D specular = spInfo.material->specular;
		// color += emissive;
		// color += ambient;
		// color += diffuse;
		// color += specular;
		return true;
	};

	processFirstIntersection( ray , BoundingBox1D( Epsilon , Infinity ) , rFilter , rKernel , spInfo , tIdx );

	// if ((ray.direction[0] - 0) < 0.1 && (ray.direction[1] + 0.707107) < 0.1 && (ray.direction[2] + 0.707107) < 0.1) {
	// 	cout << "ray in scene.todo: " << ray.direction << endl;
	// }

	return color;
}

//////////////
// Material //
//////////////
void Material::drawOpenGL( GLSLProgram * glslProgram ) const
{
	//////////////////////////////
	// Do OpenGL rendering here //
	//////////////////////////////
	WARN_ONCE( "method undefined" );

	// Sanity check to make sure that OpenGL state is good
	ASSERT_OPEN_GL_STATE();	
}

/////////////
// Texture //
/////////////
void Texture::initOpenGL( void )
{
	///////////////////////////////////
	// Do OpenGL texture set-up here //
	///////////////////////////////////
	WARN_ONCE( "method undefined" );

	// Sanity check to make sure that OpenGL state is good
	ASSERT_OPEN_GL_STATE();	
}
