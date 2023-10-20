#include <cmath>
#include <Util/exceptions.h>
#include "pointLight.h"
#include "scene.h"

using namespace Ray;
using namespace Util;
using namespace std;

////////////////
// PointLight //
////////////////

Point3D PointLight::getAmbient( Ray3D ray , const RayShapeIntersectionInfo & iInfo , const Material &material ) const
{
	////////////////////////////////////////////////////
	// Get the ambient contribution of the light here //
	////////////////////////////////////////////////////
	// WARN_ONCE( "method undefined" );
	// return Point3D();
	return _ambient + material.ambient;
}

Point3D PointLight::getDiffuse( Ray3D ray , const RayShapeIntersectionInfo &iInfo , const Material &material ) const
{
	////////////////////////////////////////////////////
	// Get the diffuse contribution of the light here //
	////////////////////////////////////////////////////
	// WARN_ONCE( "method undefined" );
	// return Point3D();
	Point3D N = iInfo.normal;
	Point3D L = _location - iInfo.position;
	Point3D K_d = material.diffuse;
	Point3D I = _diffuse; // intensity of light source after attenuation (should I be attenuating?)
	if ((ray.direction[0] - 0) < 0.1 && (ray.direction[1] + 0.707107) < 0.1 && (ray.direction[2] + 0.707107) < 0.1) {
		cout << "ray in pointLight.todo: " << ray.direction << endl;
	}
	cout << ray.direction << endl;
	return K_d * (N.dot(L)) * I;
}

Point3D PointLight::getSpecular( Ray3D ray , const RayShapeIntersectionInfo &iInfo , const Material &material ) const
{
	/////////////////////////////////////////////////////
	// Get the specular contribution of the light here //
	/////////////////////////////////////////////////////
	WARN_ONCE( "method undefined" );
	return Point3D();
}

bool PointLight::isInShadow( const RayShapeIntersectionInfo& iInfo , const Shape &shape , unsigned int tIdx ) const
{
	//////////////////////////////////////////////
	// Determine if the light is in shadow here //
	//////////////////////////////////////////////
	WARN_ONCE( "method undefined" );
	return false;
}

Point3D PointLight::transparency( const RayShapeIntersectionInfo &iInfo , const Shape &shape , Point3D cLimit , unsigned int samples , unsigned int tIdx ) const
{
	//////////////////////////////////////////////////////////
	// Compute the transparency along the path to the light //
	//////////////////////////////////////////////////////////
	WARN_ONCE( "method undefined" );
	return Point3D( 1. , 1. , 1. );
}

void PointLight::drawOpenGL( int index , GLSLProgram * glslProgram ) const
{
	//////////////////////////////
	// Do OpenGL rendering here //
	//////////////////////////////
	WARN_ONCE( "method undefined" );

	// Sanity check to make sure that OpenGL state is good
	ASSERT_OPEN_GL_STATE();	
}

