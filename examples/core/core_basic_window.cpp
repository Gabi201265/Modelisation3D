/*******************************************************************************************
*
*   raylib [core] example - Basic window
*
*   Welcome to raylib!
*
*   To test examples, just press F6 and execute raylib_compile_execute script
*   Note that compiled executable is placed in the same folder as .c file
*
*   You can find all basic examples on C:\raylib\raylib\examples folder or
*   raylib official webpage: www.raylib.com
*
*   Enjoy using raylib. :)
*
*   This example has been created using raylib 1.0 (www.raylib.com)
*   raylib is licensed under an unmodified zlib/libpng license (View raylib.h for details)
*
*   Copyright (c) 2014 Ramon Santamaria (@raysan5)
*
********************************************************************************************/

#include "raylib.h"
#include <raymath.h>
#include "rlgl.h"
#include <math.h>
#include <float.h>
#include <vector>

#if defined(PLATFORM_DESKTOP)
#define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
#define GLSL_VERSION            100
#endif

#define EPSILON 1.e-6f
#define GRAVITY Vector3 {0,-9.82f,0}

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}


struct Cylindrical { //Définition des coordonnées cylindrique
	float rho;
	float theta;
	float y;

	inline Cylindrical operator+(Cylindrical a) {
		return { a.rho + rho,a.theta + theta,a.y + y };
	}
};
struct Spherical { //Définition des coordonnées sphérique
	float rho;
	float theta;
	float phi;

	inline Spherical operator+(Spherical a) {
		return { a.rho + rho,a.theta + theta,a.phi + phi };
	}

};
struct Segment {
	Vector3 pt1;
	Vector3 pt2;
};

struct Plane {
	Vector3 normal;
	float d;//Distance signée
};

struct Sphere {
	Vector3 center;//centre
	float r;//Rayon
};

struct Cylinder {
	Vector3 pt1;
	Vector3 pt2;
	float radius;
};

struct Referential {//Structure de chaque référentiel
	Vector3 origin;
	Vector3 i, j, k;
	Quaternion qua;

};

struct Disc {
	Referential ref;
	float radius;
};

struct Quad {
	Referential ref;
	//Demi size
	Vector2 demSize;
};

struct Capsule {
	Vector3 pt1;
	Vector3 pt2;
	float radius;
	Quaternion qua;
};

struct Box{
	Referential ref;
	//Demi size
	Vector3 demSize;
};

struct Box_arrondie{
	float rayon;
	Box box;
	Capsule listCapsule[12];
	Quad listQuad[6];
};

struct Ball {
	float rayon;
	Vector3 position;
	Vector3 mvmt;
};

inline void vertex(Vector3 v) {
	rlVertex3f(v.x, v.y, v.z);
}

Cylindrical CartesianToCylindrical(Vector3 cart)//Conversion voir feuille
{
	Cylindrical cyl;
	cyl.rho = sqrtf(cart.x * cart.x + cart.z * cart.z);
	cyl.y = cart.y;

	if (cyl.rho < EPSILON)cyl.theta = 0;
	else
	{
		cyl.theta = atan2f(cart.x, cart.z);
		if (cyl.theta < 0)cyl.theta += PI * 2;
	}
	return cyl;
}

Vector3 CylindricalToCartesian(Cylindrical cyl)//Conversion voir feuille
{
	return Vector3{ cyl.rho * sinf(cyl.theta),cyl.y,cyl.rho * cosf(cyl.theta) };
}

Vector3 SphericalToCartesian(Spherical sph)//Conversion voir feuille
{
	return Vector3{ sph.rho * sinf(sph.phi) * sinf(sph.theta),
	sph.rho * cosf(sph.phi),
	sph.rho * sinf(sph.phi) * cosf(sph.theta) };
}

Vector3 GlobalToLocalPos(Vector3 globalPos, Referential localRef) {
	Vector3 oPrimeP = Vector3Subtract(globalPos, localRef.origin);
	Vector3 vec = { Vector3DotProduct(oPrimeP, localRef.i) , Vector3DotProduct(oPrimeP, localRef.j) , Vector3DotProduct(oPrimeP, localRef.k) };
	return vec;
}

Vector3 GlobalToLocalVect(Vector3 globalVect, Referential localRef) {
	Vector3 vec = { Vector3DotProduct(globalVect, localRef.i) , Vector3DotProduct(globalVect, localRef.j) , Vector3DotProduct(globalVect, localRef.k) };
	return vec;
}

Vector3 LocalToGlobalPos(Vector3 localPos, Referential localRef) {
	Vector3 oPrimeP = Vector3Add(Vector3Add(Vector3Scale(localRef.i, localPos.x), Vector3Scale(localRef.j, localPos.y)), Vector3Scale(localRef.k, localPos.z));
	return Vector3Add(localRef.origin, oPrimeP);
}

Vector3 LocalToGlobalVect(Vector3 localVect, Referential localRef) {
	Vector3 vecx = Vector3Scale(localRef.i, localVect.x);
	Vector3 vecy = Vector3Scale(localRef.j, localVect.y);
	Vector3 vecz = Vector3Scale(localRef.k, localVect.z);
	return Vector3Add(Vector3Add(vecx, vecy), vecz);
}


void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
	static Spherical sphPos = { 10,PI / 4.f,PI / 4.f };
	static Spherical sphSpeed = { 10,.4f,.4f };
	float rhoMin = 4;
	float rhoMax = 40;

	static Vector2 prevMousePos = { 0,0 };
	Vector2 mousePos = GetMousePosition();
	Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
	prevMousePos = mousePos;

	Spherical sphDelta = { -GetMouseWheelMove() * sphSpeed.rho * deltaTime, //Calculer en spérique une PETITE variation (delta) pour ro (-le mouvement de la souris * la vitesse * le temps
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.x * sphSpeed.theta * deltaTime : 0, //Pour theta seulement si le bouton droit de la souris est cliqué (mouvement horizontal)
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.y * sphSpeed.phi * deltaTime : 0 }; //Pour phi seulement si le bouton droit de la souris est cliqué (mouvement vertical)

	Spherical newSphPos = sphPos + sphDelta; //La nouvelle position sphérique est donc la position initiale et la variation calculer ci-avant
	newSphPos = { Clamp(newSphPos.rho,rhoMin,rhoMax), //Clamp dans un intervalle de donnée numérique. Si dans l'intervalle alors elle conserve sa valeur, sinon elle prend la valeur de la borne la plus proche. On le fait pour ro
	newSphPos.theta, //Pas besoin de clamp car on a dit que theta n'avait pas de limite
	Clamp(newSphPos.phi,PI / 100.f,.99f * PI) }; //On le clamp car limité (phi). Je n'autorise pas la caméra a avoir une vue plongeante complètement vers le haut ou complètement vers le bas. On perd la caméra, elle ne sait plus comment s'orienter par rapport à l'axe "guimboloque"

	sphPos = newSphPos; //La position courante est donc la nouvelle position de la caméra

	camera->position = SphericalToCartesian(sphPos); //On fixe la position de la caméra en convertissant les coordonnées sphériques en cartésien
}

Quad* ListQuadrilatere(Box_arrondie roundedbox) {
	//Vector3RotateByQuaternion(roundedbox.box.ref.i,roundedbox.box.ref.qua)
	Referential ref;
	Quaternion rot;
	roundedbox.listQuad[0] =  { {Vector3Add(roundedbox.box.ref.origin, LocalToGlobalVect({ 0, roundedbox.box.demSize.y + roundedbox.rayon, 0 }, roundedbox.box.ref)), roundedbox.box.ref.i,roundedbox.box.ref.j,roundedbox.box.ref.k, roundedbox.box.ref.qua}, { roundedbox.box.demSize.x, roundedbox.box.demSize.z } }; //Au dessus

	rot = QuaternionFromAxisAngle(roundedbox.box.ref.k, PI);
	ref = { Vector3Add(roundedbox.box.ref.origin,LocalToGlobalVect({ 0, -(roundedbox.box.demSize.y + roundedbox.rayon),0},roundedbox.box.ref)),Vector3RotateByQuaternion(roundedbox.box.ref.i,rot),roundedbox.box.ref.j,roundedbox.box.ref.k,QuaternionMultiply(rot,roundedbox.box.ref.qua)};
	roundedbox.listQuad[1] = {ref, { roundedbox.box.demSize.x, roundedbox.box.demSize.z }}; //En dessous
	
	rot = QuaternionFromAxisAngle(roundedbox.box.ref.k, -PI / 2);
	ref = { Vector3Add(roundedbox.box.ref.origin,LocalToGlobalVect({roundedbox.box.demSize.x + roundedbox.rayon, 0, 0},roundedbox.box.ref)),Vector3RotateByQuaternion(roundedbox.box.ref.i,rot),Vector3RotateByQuaternion(roundedbox.box.ref.j,rot),Vector3RotateByQuaternion(roundedbox.box.ref.k,rot),QuaternionMultiply(rot,roundedbox.box.ref.qua) };
	roundedbox.listQuad[2] = {ref, {roundedbox.box.demSize.y, roundedbox.box.demSize.z}}; //A droite
	
	rot = QuaternionFromAxisAngle(roundedbox.box.ref.k, PI / 2);
	ref = { Vector3Add(roundedbox.box.ref.origin,LocalToGlobalVect({-(roundedbox.box.demSize.x + roundedbox.rayon), 0, 0},roundedbox.box.ref)),Vector3RotateByQuaternion(roundedbox.box.ref.i,rot),Vector3RotateByQuaternion(roundedbox.box.ref.j,rot),Vector3RotateByQuaternion(roundedbox.box.ref.k,rot),QuaternionMultiply(rot,roundedbox.box.ref.qua) };
	roundedbox.listQuad[3] = {ref, {roundedbox.box.demSize.y, roundedbox.box.demSize.z}}; //A gauche
	
	rot = QuaternionFromAxisAngle(roundedbox.box.ref.i, PI / 2);
	ref = { Vector3Add(roundedbox.box.ref.origin,LocalToGlobalVect({0, 0, roundedbox.box.demSize.z + roundedbox.rayon},roundedbox.box.ref)),Vector3RotateByQuaternion(roundedbox.box.ref.i,rot),Vector3RotateByQuaternion(roundedbox.box.ref.j,rot),Vector3RotateByQuaternion(roundedbox.box.ref.k,rot),QuaternionMultiply(rot,roundedbox.box.ref.qua) };
	roundedbox.listQuad[4] = {ref, {roundedbox.box.demSize.x, roundedbox.box.demSize.y}}; //De face
	
	rot = QuaternionFromAxisAngle(roundedbox.box.ref.i, -PI / 2);
	ref= { Vector3Add(roundedbox.box.ref.origin,LocalToGlobalVect({0, 0, -(roundedbox.box.demSize.z + roundedbox.rayon)},roundedbox.box.ref)),Vector3RotateByQuaternion(roundedbox.box.ref.i,rot),Vector3RotateByQuaternion(roundedbox.box.ref.j,rot),Vector3RotateByQuaternion(roundedbox.box.ref.k,rot),QuaternionMultiply(rot,roundedbox.box.ref.qua) };
	roundedbox.listQuad[5] = {ref, { roundedbox.box.demSize.x, roundedbox.box.demSize.y } }; //De dos
	return roundedbox.listQuad;
}

Capsule* ListCapsule(Box_arrondie boxrounded) {

	Vector3 bottomFrontLeft = Vector3Add(Vector3Add(Vector3Add(boxrounded.box.ref.origin, Vector3Scale(boxrounded.box.ref.i, -boxrounded.box.demSize.x)), Vector3Scale(boxrounded.box.ref.j, -boxrounded.box.demSize.y)), Vector3Scale(boxrounded.box.ref.k, boxrounded.box.demSize.z));
	Vector3 bottomFrontRight = Vector3Add(Vector3Add(Vector3Add(boxrounded.box.ref.origin, Vector3Scale(boxrounded.box.ref.i, boxrounded.box.demSize.x)), Vector3Scale(boxrounded.box.ref.j, -boxrounded.box.demSize.y)), Vector3Scale(boxrounded.box.ref.k, boxrounded.box.demSize.z));
	Vector3 bottomBackLeft = Vector3Add(Vector3Add(Vector3Add(boxrounded.box.ref.origin, Vector3Scale(boxrounded.box.ref.i, -boxrounded.box.demSize.x)), Vector3Scale(boxrounded.box.ref.j, -boxrounded.box.demSize.y)), Vector3Scale(boxrounded.box.ref.k, -boxrounded.box.demSize.z));
	Vector3 bottomBackRight = Vector3Add(Vector3Add(Vector3Add(boxrounded.box.ref.origin, Vector3Scale(boxrounded.box.ref.i, boxrounded.box.demSize.x)), Vector3Scale(boxrounded.box.ref.j, -boxrounded.box.demSize.y)), Vector3Scale(boxrounded.box.ref.k, -boxrounded.box.demSize.z));
	Vector3 topFrontLeft = Vector3Add(Vector3Add(Vector3Add(boxrounded.box.ref.origin, Vector3Scale(boxrounded.box.ref.i, -boxrounded.box.demSize.x)), Vector3Scale(boxrounded.box.ref.j, boxrounded.box.demSize.y)), Vector3Scale(boxrounded.box.ref.k, boxrounded.box.demSize.z));
	Vector3 topFrontRight = Vector3Add(Vector3Add(Vector3Add(boxrounded.box.ref.origin, Vector3Scale(boxrounded.box.ref.i, boxrounded.box.demSize.x)), Vector3Scale(boxrounded.box.ref.j, boxrounded.box.demSize.y)), Vector3Scale(boxrounded.box.ref.k, boxrounded.box.demSize.z));
	Vector3 topBackLeft = Vector3Add(Vector3Add(Vector3Add(boxrounded.box.ref.origin, Vector3Scale(boxrounded.box.ref.i, -boxrounded.box.demSize.x)), Vector3Scale(boxrounded.box.ref.j, boxrounded.box.demSize.y)), Vector3Scale(boxrounded.box.ref.k, -boxrounded.box.demSize.z));
	Vector3 topBackRight = Vector3Add(Vector3Add(Vector3Add(boxrounded.box.ref.origin, Vector3Scale(boxrounded.box.ref.i, boxrounded.box.demSize.x)), Vector3Scale(boxrounded.box.ref.j, boxrounded.box.demSize.y)), Vector3Scale(boxrounded.box.ref.k, -boxrounded.box.demSize.z));

	Capsule caps1 = { bottomFrontLeft,bottomFrontRight,boxrounded.rayon, QuaternionMultiply(boxrounded.box.ref.qua,QuaternionFromAxisAngle({ 0,0,1 }, -PI / 2)) };
	boxrounded.listCapsule[0] = caps1;
	Capsule caps2 = { topFrontLeft,topFrontRight,boxrounded.rayon,QuaternionMultiply(boxrounded.box.ref.qua,QuaternionFromAxisAngle({ 0,0,1 }, -PI / 2)) };
	boxrounded.listCapsule[1] = caps2;
	Capsule caps3 = { topFrontLeft,bottomFrontLeft,boxrounded.rayon,QuaternionMultiply(boxrounded.box.ref.qua,QuaternionFromAxisAngle({ 1,0,0 }, -PI)) };
	boxrounded.listCapsule[2] = caps3;
	Capsule caps4 = { topFrontRight,bottomFrontRight,boxrounded.rayon,QuaternionMultiply(boxrounded.box.ref.qua,QuaternionFromAxisAngle({ 1,0,0 }, -PI)) };
	boxrounded.listCapsule[3] = caps4;

	//Back face
	Capsule caps5 = { topBackRight,topBackLeft,boxrounded.rayon, QuaternionMultiply(boxrounded.box.ref.qua,QuaternionFromAxisAngle({ 0,0,1 }, PI/2)) };
	boxrounded.listCapsule[4] = caps5;
	Capsule caps6 = { bottomBackLeft,bottomBackRight,boxrounded.rayon,QuaternionMultiply(boxrounded.box.ref.qua,QuaternionFromAxisAngle({ 0,0,1 }, -PI / 2)) };
	boxrounded.listCapsule[5] = caps6;
	Capsule caps7 = { topBackLeft,bottomBackLeft,boxrounded.rayon,QuaternionMultiply(boxrounded.box.ref.qua,QuaternionFromAxisAngle({ 0,0,1 }, PI)) };
	boxrounded.listCapsule[6] = caps7;
	Capsule caps8 = { topBackRight,bottomBackRight,boxrounded.rayon,QuaternionMultiply(boxrounded.box.ref.qua,QuaternionFromAxisAngle({ 0,0,1 }, PI)) };
	boxrounded.listCapsule[7] = caps8;

	//Right face
	Capsule caps9 = { topFrontRight,topBackRight,boxrounded.rayon,QuaternionMultiply(QuaternionFromAxisAngle({ 1,0,0 }, -PI/2),boxrounded.box.ref.qua) };//TopLeft et TopRight
	boxrounded.listCapsule[8] = caps9;
	Capsule caps10 = { bottomFrontRight ,bottomBackRight,boxrounded.rayon,QuaternionMultiply(QuaternionFromAxisAngle({ 1,0,0 }, -PI / 2),boxrounded.box.ref.qua) };//BttomLeft et BottomRight
	boxrounded.listCapsule[9] = caps10;

	//Left face
	Capsule caps11 = { topFrontLeft,topBackLeft,boxrounded.rayon ,QuaternionMultiply(QuaternionFromAxisAngle({ 1,0,0 }, -PI / 2),boxrounded.box.ref.qua) };//TopLeft et TopRight
	boxrounded.listCapsule[10] = caps11;
	Capsule caps12 = { bottomFrontLeft,bottomBackLeft,boxrounded.rayon,QuaternionMultiply(QuaternionFromAxisAngle({ 1,0,0 }, -PI / 2),boxrounded.box.ref.qua) };//BottomRight et TopRight
	boxrounded.listCapsule[11] = caps12;


	return boxrounded.listCapsule;
}

//Méthodes de dessins et d'affichages

// Draw sphere with extended parameters
void MyDrawSphereEx2(Quaternion q, Vector3 centerPos, Vector3 radius, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,1,0 });

	int numVertex = nSegmentsTheta * nSegmentsPhi * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(centerPos.x, centerPos.y, centerPos.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlScalef(radius.x, radius.y, radius.z);


	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaPhi = PI / nSegmentsPhi;
	float deltaTheta = 2 * PI / nSegmentsTheta;

	float phi = 0;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = 0;
		Vector3 tmpBottomLeft = SphericalToCartesian(Spherical{ 1,theta,phi + deltaPhi });

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = tmpBottomLeft;
			Vector3 topRight = vertexBufferTheta[j + 1];
			Vector3 bottomRight = SphericalToCartesian(Spherical{ 1,theta + deltaTheta,phi + deltaPhi });


			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);
			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			theta += deltaTheta;

			vertexBufferTheta[j] = tmpBottomLeft;
			tmpBottomLeft = bottomRight;
		}
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}
	rlEnd();
	rlPopMatrix();
}

void MyDrawSphereWiresEx2(Quaternion q, Vector3 centerPos, Vector3 radius, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,1,0 });

	int numVertex = nSegmentsTheta * nSegmentsPhi * 4;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(centerPos.x, centerPos.y, centerPos.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);


	rlScalef(radius.x, radius.y, radius.z);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaPhi = PI / nSegmentsPhi;
	float deltaTheta = 2 * PI / nSegmentsTheta;

	float phi = 0;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = 0;

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = SphericalToCartesian(Spherical{ 1,theta,phi + deltaPhi });
			Vector3 topRight = vertexBufferTheta[j + 1];

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);

			theta += deltaTheta;

			vertexBufferTheta[j] = bottomLeft;
		}
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}
	rlEnd();
	rlPopMatrix();
}


void MyDrawSphereTrianglesAndWiresEx(Vector3 centerPos, float radius, int nSegmentsTheta, int nSegmentsPhi, Color trianglesColor, Color wiresColor)
{
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,radius,0 });

	int numVertex = nSegmentsTheta * nSegmentsPhi * 10;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(centerPos.x, centerPos.y, centerPos.z);
	rlScalef(radius, radius, radius);


	float deltaPhi = PI / nSegmentsPhi;
	float deltaTheta = 2 * PI / nSegmentsTheta;

	float phi = 0;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = 0;
		Vector3 tmpBottomLeft = SphericalToCartesian(Spherical{ radius,theta,phi + deltaPhi });

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = tmpBottomLeft;
			Vector3 topRight = vertexBufferTheta[j + 1];
			Vector3 bottomRight = SphericalToCartesian(Spherical{ radius,theta + deltaTheta,phi + deltaPhi });

			rlBegin(RL_TRIANGLES);
			rlColor4ub(trianglesColor.r, trianglesColor.g, trianglesColor.b, trianglesColor.a);

			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);
			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);
			rlEnd();


			rlBegin(RL_LINES);
			rlColor4ub(wiresColor.r, wiresColor.g, wiresColor.b, wiresColor.a);
			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);
			rlEnd();


			theta += deltaTheta;

			vertexBufferTheta[j] = tmpBottomLeft;
			tmpBottomLeft = bottomRight;
		}
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}

	rlPopMatrix();
}


void MyDrawQuad(Quaternion q, Vector3 centre, Vector2 size, Color color) {
	int numVertex = 6;
	if (rlCheckBufferLimit(numVertex))
		rlglDraw();

	rlPushMatrix();
	rlTranslatef(centre.x, centre.y, centre.z);
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(size.x, 0, size.y);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	rlVertex3f(-1, 0, -1);
	rlVertex3f(-1, 0, 1);
	rlVertex3f(1, 0, -1);

	rlVertex3f(1, 0, -1);
	rlVertex3f(-1, 0, 1);
	rlVertex3f(1, 0, 1);

	rlEnd();
	rlPopMatrix();
}


void MyDrawQuadWire(Quaternion q, Vector3 center, Vector2 size, Color color) {

	rlglDraw();
	rlPushMatrix();
	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	Vector3 topLeft = { center.x - (size.x / 2) , center.y, center.z - (size.y / 2) };
	Vector3 bottomLeft = { center.x - (size.x / 2), center.y , center.z + (size.y / 2) };
	Vector3 topRight = { center.x + (size.x / 2), center.y , center.z - (size.y / 2) };
	Vector3 bottomRight = { center.x + (size.x / 2), center.y, center.z + (size.y / 2) };

	DrawLine3D(topLeft, bottomRight, color);//Méthode de raylib
	DrawLine3D(topRight, bottomLeft, color);
	DrawLine3D(topLeft, topRight, color);
	DrawLine3D(bottomLeft, bottomRight, color);
	DrawLine3D(topLeft, bottomLeft, color);
	DrawLine3D(topRight, bottomRight, color);
	rlEnd();
	rlPopMatrix();
}
void MyDrawDisk(Quaternion q, Vector3 center, float radius, int nSegmentsTheta, Color color) {
	if (nSegmentsTheta < 3) return;

	int numVertex = nSegmentsTheta * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	rlPushMatrix();

	rlTranslatef(center.x, center.y, center.z);//bouger matrice entière, tous les vecteurs, la couleur...comme un seul bloc

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	float deltaTheta = 2 * PI / nSegmentsTheta;
	for (float i = 0; i < 2 * PI; i += deltaTheta)
	{

		Vector3 pt1 = CylindricalToCartesian(Cylindrical{ radius,i,0 });
		Vector3 pt2 = CylindricalToCartesian(Cylindrical{ radius,i + deltaTheta,0 });

		rlVertex3f(pt1.x, pt1.y, pt1.z);
		//rlVertex3f(center.x, center.y, center.z);
		rlVertex3f(pt2.x, pt2.y, pt2.z);
		rlVertex3f(0, 0, 0);
	}
	rlEnd();
	rlPopMatrix();

}
void MyDrawCylinder(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color) {
	// le paramètre booléen drawCaps sert à déclencher l'affichage des extrémités discoïdes du cylindre
	// dans le cas d'un cylindre infini, les extrémités discoïdes ne sont pas affichées, donc drawCaps = false
	// dans le cas d'un cylindre fini,  les extrémités discoïdes sont affichées, donc drawCaps = true
	if (nSegmentsTheta < 3) return;

	int numVertex = nSegmentsTheta * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	rlPushMatrix();

	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);//bouger matrice entière, tous les vecteurs, la couleur...comme un seul bloc

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	float height = Vector3Distance(cyl.pt1, cyl.pt2);;

	rlScalef(cyl.radius, height, cyl.radius);    //Défini le rayon du cylindre! (sur X et Z), defini également la longueur du cylindre (sur Y)

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	float deltaTheta = 2 * PI / nSegmentsTheta;
	if (drawCaps == true) {

		MyDrawDisk(QuaternionFromAxisAngle({ 1,0,0 }, PI), { 0,0,0 }, 1, nSegmentsTheta, color);
		MyDrawDisk(QuaternionFromAxisAngle({ 1,0,0 }, 0), { 0,1,0 }, 1, nSegmentsTheta, color);
	}
	for (float i = 0; i < 2 * PI; i += deltaTheta)
	{
		Vector3 pt1 = CylindricalToCartesian(Cylindrical{ 1,i,0 });//Bottom Left
		Vector3 pt2 = CylindricalToCartesian(Cylindrical{ 1,i + deltaTheta,0 });//Bottom Right
		Vector3 pt3 = CylindricalToCartesian(Cylindrical{ 1,i + deltaTheta,1 });//Top Right
		rlVertex3f(pt1.x, pt1.y, pt1.z);
		rlVertex3f(pt2.x, pt2.y, pt2.z);
		rlVertex3f(pt3.x, pt3.y, pt3.z);

		Vector3 pt4 = CylindricalToCartesian(Cylindrical{ 1,i,1 });//Top Left
		Vector3 pt5 = CylindricalToCartesian(Cylindrical{ 1,i,0 });//Bottom Left
		Vector3 pt6 = CylindricalToCartesian(Cylindrical{ 1,i + deltaTheta,1 });//Top Right
		rlVertex3f(pt4.x, pt4.y, pt4.z);
		rlVertex3f(pt5.x, pt5.y, pt5.z);
		rlVertex3f(pt6.x, pt6.y, pt6.z);
	}
	rlEnd();
	rlPopMatrix();


}

void MyDrawCaps(Quaternion q, Capsule caps, int nSegmentsTheta, Color color) {
	//On recopie MyDrawCylinder en changeant les extrémités
	if (nSegmentsTheta < 3) return;
	int numVertex = nSegmentsTheta * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	rlTranslatef(caps.pt1.x, caps.pt1.y, caps.pt1.z);//bouger matrice entière, tous les vecteurs, la couleur...comme un seul bloc

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	float height = Vector3Distance(caps.pt1, caps.pt2);

	rlScalef(caps.radius, height, caps.radius);									//Défini le rayon du cylindre! (sur X et Z), defini également la longueur du cylindre (sur Y)
																				//DrawSphere(caps.pt1, caps.radius, color);

	MyDrawSphereEx2(QuaternionFromAxisAngle({ 1,0,0 }, PI), { 0,0,0 }, Vector3{ 1, caps.radius / height, 1 }, nSegmentsTheta, 20, color);
	MyDrawSphereWiresEx2(QuaternionFromAxisAngle({ 1,0,0 }, PI), { 0,0,0 }, Vector3{ 1, caps.radius / height, 1 }, nSegmentsTheta, 20, WHITE);
	MyDrawSphereEx2(QuaternionFromAxisAngle({ 1,0,0 }, 0), { 0,1,0 }, Vector3{ 1, caps.radius / height, 1 }, nSegmentsTheta, 20, color);
	MyDrawSphereWiresEx2(QuaternionFromAxisAngle({ 1,0,0 }, 0), { 0,1,0 }, Vector3{ 1, caps.radius / height, 1 }, nSegmentsTheta, 20, WHITE);

	rlBegin(RL_TRIANGLES);

	rlColor4ub(color.r, color.g, color.b, color.a);
	float deltaTheta = 2 * PI / nSegmentsTheta;
	for (float i = 0; i < 2 * PI; i += deltaTheta)
	{
		Vector3 pt1 = CylindricalToCartesian(Cylindrical{ 1,i,0 });//Bottom Left
		Vector3 pt2 = CylindricalToCartesian(Cylindrical{ 1,i + deltaTheta,0 });//Bottom Right
		Vector3 pt3 = CylindricalToCartesian(Cylindrical{ 1,i + deltaTheta, 1 });//Top Right
		rlVertex3f(pt1.x, pt1.y, pt1.z);
		rlVertex3f(pt2.x, pt2.y, pt2.z);
		rlVertex3f(pt3.x, pt3.y, pt3.z);

		Vector3 pt4 = CylindricalToCartesian(Cylindrical{ 1,i, 1 });//Top Left
		Vector3 pt5 = CylindricalToCartesian(Cylindrical{ 1,i,0 });//Bottom Left
		Vector3 pt6 = CylindricalToCartesian(Cylindrical{ 1,i + deltaTheta, 1 });//Top Right
		rlVertex3f(pt4.x, pt4.y, pt4.z);
		rlVertex3f(pt5.x, pt5.y, pt5.z);
		rlVertex3f(pt6.x, pt6.y, pt6.z);
	}

	rlEnd();
	rlPopMatrix();
}
void MyDrawWiresBox(Quaternion q, Box box, Color color) {

	Vector3 bottomFrontLeft = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, -box.demSize.x)), Vector3Scale(box.ref.j, -box.demSize.y)), Vector3Scale(box.ref.k, box.demSize.z));
	Vector3 bottomFrontRight = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, box.demSize.x)), Vector3Scale(box.ref.j, -box.demSize.y)), Vector3Scale(box.ref.k, box.demSize.z));
	Vector3 bottomBackLeft = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, -box.demSize.x)), Vector3Scale(box.ref.j, -box.demSize.y)), Vector3Scale(box.ref.k, -box.demSize.z));
	Vector3 bottomBackRight = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, box.demSize.x)), Vector3Scale(box.ref.j, -box.demSize.y)), Vector3Scale(box.ref.k, -box.demSize.z));
	Vector3 topFrontLeft = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, -box.demSize.x)), Vector3Scale(box.ref.j, box.demSize.y)), Vector3Scale(box.ref.k, box.demSize.z));
	Vector3 topFrontRight = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, box.demSize.x)), Vector3Scale(box.ref.j, box.demSize.y)), Vector3Scale(box.ref.k, box.demSize.z));
	Vector3 topBackLeft = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, -box.demSize.x)), Vector3Scale(box.ref.j, box.demSize.y)), Vector3Scale(box.ref.k, -box.demSize.z));
	Vector3 topBackRight = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, box.demSize.x)), Vector3Scale(box.ref.j, box.demSize.y)), Vector3Scale(box.ref.k, -box.demSize.z));
	if (rlCheckBufferLimit(36)) rlglDraw();

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> rotate -> translate)

	//rlTranslatef(position.x, position.y, position.z);

	//Scalef(1.0f, 1.0f, 1.0f);

		rlBegin(RL_LINES);
		rlColor4ub(color.r, color.g, color.b, color.a);
		//Face de devant
		vertex(bottomFrontLeft);
		vertex(bottomFrontRight);

		vertex(topFrontLeft);
		vertex(bottomFrontRight);

		vertex(topFrontRight);
		vertex(topFrontLeft);

		//Face du dessus
		vertex(topBackLeft);
		vertex(topFrontLeft);

		vertex(topBackRight);
		vertex(topFrontLeft);

		vertex(topFrontRight);
		vertex(topBackRight);

		//face coté droit
		vertex(topBackRight);
		vertex(topFrontRight);

		vertex(bottomBackRight);
		vertex(topFrontRight);

		vertex(bottomFrontRight);
		vertex(bottomBackRight);

		//Face coté gauche
		vertex(topFrontLeft);
		vertex(topBackLeft);

		vertex(bottomBackLeft);
		vertex(topFrontLeft);

		vertex(bottomBackLeft);
		vertex(bottomFrontLeft);


		//Face du dessous
		vertex(bottomFrontLeft);
		vertex(bottomBackLeft);

		vertex(bottomBackRight);
		vertex(bottomFrontLeft);

		vertex(bottomBackRight);
		vertex(bottomFrontRight);


		//Face de derriere
		vertex(bottomBackLeft);
		vertex(topBackRight);

		vertex(bottomBackRight);
		vertex(bottomBackLeft);

		vertex(topBackLeft);
		vertex(topBackRight);

		rlEnd();
	rlPopMatrix();

}
void MyDrawBox(Quaternion q, Box box, Color color) {

	Vector3 bottomFrontLeft = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, -box.demSize.x)), Vector3Scale(box.ref.j, -box.demSize.y)), Vector3Scale(box.ref.k, box.demSize.z));
	Vector3 bottomFrontRight = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, box.demSize.x)), Vector3Scale(box.ref.j, -box.demSize.y)), Vector3Scale(box.ref.k, box.demSize.z));
	Vector3 bottomBackLeft = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, -box.demSize.x)), Vector3Scale(box.ref.j, -box.demSize.y)), Vector3Scale(box.ref.k, -box.demSize.z));
	Vector3 bottomBackRight = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, box.demSize.x)), Vector3Scale(box.ref.j, -box.demSize.y)), Vector3Scale(box.ref.k, -box.demSize.z));
	Vector3 topFrontLeft = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, -box.demSize.x)), Vector3Scale(box.ref.j, box.demSize.y)), Vector3Scale(box.ref.k, box.demSize.z));
	Vector3 topFrontRight = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, box.demSize.x)), Vector3Scale(box.ref.j, box.demSize.y)), Vector3Scale(box.ref.k, box.demSize.z));
	Vector3 topBackLeft = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, -box.demSize.x)), Vector3Scale(box.ref.j, box.demSize.y)), Vector3Scale(box.ref.k, -box.demSize.z));
	Vector3 topBackRight = Vector3Add(Vector3Add(Vector3Add(box.ref.origin, Vector3Scale(box.ref.i, box.demSize.x)), Vector3Scale(box.ref.j, box.demSize.y)), Vector3Scale(box.ref.k, -box.demSize.z));
	
	if (rlCheckBufferLimit(36)) rlglDraw();
	
	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> rotate -> translate)
	
	//rlTranslatef(position.x, position.y, position.z);

	//Scalef(1.0f, 1.0f, 1.0f);

		rlBegin(RL_TRIANGLES);
		rlColor4ub(color.r, color.g, color.b, color.a);

		//Face de devant
		vertex(bottomFrontLeft);
		vertex(bottomFrontRight);

		vertex(topFrontLeft);
		vertex(bottomFrontRight);

		vertex(topFrontRight);
		vertex(topFrontLeft);

		//Face du dessus
		vertex(topBackLeft);
		vertex(topFrontLeft);

		vertex(topBackRight);
		vertex(topFrontLeft);

		vertex(topFrontRight);
		vertex(topBackRight);

		//face coté droit
		vertex(topBackRight);
		vertex(topFrontRight);

		vertex(bottomBackRight);
		vertex(topFrontRight);

		vertex(bottomFrontRight);
		vertex(bottomBackRight);

		//Face coté gauche
		vertex(topFrontLeft);
		vertex(topBackLeft);

		vertex(bottomBackLeft);
		vertex(topFrontLeft);

		vertex(bottomBackLeft);
		vertex(bottomFrontLeft);
		
		
		//Face du dessous
		vertex(bottomFrontLeft);
		vertex(bottomBackLeft);

		vertex(bottomBackRight);
		vertex(bottomFrontLeft);

		vertex(bottomBackRight);
		vertex(bottomFrontRight);
		

		//Face de derriere
		vertex(bottomBackLeft);
		vertex(topBackRight);

		vertex(bottomBackRight);
		vertex(bottomBackLeft);

		vertex(topBackLeft);
		vertex(topBackRight);

		rlEnd();
	rlPopMatrix();
}

void MyDrawRoundedBox(Quaternion q, Box_arrondie boxrounded, Color color) {
	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;

	//Mettre le quaternion a 0 sinon il s'applique deux fois
	Quad* tab = ListQuadrilatere(boxrounded);
	MyDrawQuad(tab[0].ref.qua, tab[0].ref.origin, tab[0].demSize, RED);
	MyDrawQuad(tab[1].ref.qua, tab[1].ref.origin, tab[1].demSize, RED);
	MyDrawQuad(tab[2].ref.qua, tab[2].ref.origin, tab[2].demSize, RED);
	MyDrawQuad(tab[3].ref.qua, tab[3].ref.origin, tab[3].demSize, RED);
	MyDrawQuad(tab[4].ref.qua, tab[4].ref.origin, tab[4].demSize, RED);
	MyDrawQuad(tab[5].ref.qua, tab[5].ref.origin, tab[5].demSize, RED);

	if (rlCheckBufferLimit(36)) rlglDraw();
	rlPushMatrix();
		rlTranslatef(boxrounded.box.ref.origin.x, boxrounded.box.ref.origin.y, boxrounded.box.ref.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	
	//Front face
	Capsule caps1 = {{x - boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z + boxrounded.box.demSize.z},{x - boxrounded.box.demSize.x, y + boxrounded.box.demSize.y , z + boxrounded.box.demSize.z},boxrounded.rayon,QuaternionFromAxisAngle({1,0,0}, 0)};//BottomLeft et TopLeft
	MyDrawCaps(caps1.qua, caps1, 60, BLACK);
	Capsule caps2 = { {x - boxrounded.box.demSize.x , y + boxrounded.box.demSize.y,z + boxrounded.box.demSize.z },{x + boxrounded.box.demSize.x , y + boxrounded.box.demSize.y, z + boxrounded.box.demSize.z },boxrounded.rayon,QuaternionFromAxisAngle({ 0,0,1 }, -PI / 2) };//TopLeft et TopRight
	MyDrawCaps(caps2.qua, caps2, 60, BLACK);
	Capsule caps3 = { {x - boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z + boxrounded.box.demSize.z },{x + boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z + boxrounded.box.demSize.z },boxrounded.rayon,QuaternionFromAxisAngle({ 0,0,1 }, -PI / 2) };//BottomLeft et BottomRight
	MyDrawCaps(caps3.qua , caps3, 60, BLACK);
	Capsule caps4 = { {x + boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z + boxrounded.box.demSize.z },{x + boxrounded.box.demSize.x , y + boxrounded.box.demSize.y , z + boxrounded.box.demSize.z },boxrounded.rayon, QuaternionFromAxisAngle({ 1,0,0 }, 0) };//BottomRight et TopRight
	MyDrawCaps(caps4.qua, caps4, 60, BLACK);

	//Back face
	Capsule caps5 = { {x - boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },{x - boxrounded.box.demSize.x , y + boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },boxrounded.rayon, QuaternionFromAxisAngle({ 1,0,0 }, 0) };//BottomLeft et TopLeft
	MyDrawCaps(caps5.qua, caps5, 60, BLACK);
	Capsule caps6 = { {x - boxrounded.box.demSize.x , y + boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },{x + boxrounded.box.demSize.x , y + boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },boxrounded.rayon,QuaternionFromAxisAngle({ 0,0,1 }, -PI / 2) };//TopLeft et TopRight
	MyDrawCaps(caps6.qua, caps6, 60, BLACK);
	Capsule caps7 = { {x - boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },{x + boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },boxrounded.rayon,QuaternionFromAxisAngle({ 0,0,1 }, -PI / 2) };//BottomLeft et BottomRight
	MyDrawCaps(caps7.qua, caps7, 60, BLACK);
	Capsule caps8 = { {x + boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },{x + boxrounded.box.demSize.x , y + boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },boxrounded.rayon,QuaternionFromAxisAngle({ 1,0,0 }, 0) };//BottomRight et TopRight
	MyDrawCaps(caps8.qua, caps8, 60, BLACK);

	//Right face
	Capsule caps9 = { {x - boxrounded.box.demSize.x , y + boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },{x - boxrounded.box.demSize.x , y + boxrounded.box.demSize.y , z + boxrounded.box.demSize.z },boxrounded.rayon,QuaternionFromAxisAngle({ 1,0,0 }, PI / 2) };//TopLeft et TopRight
	MyDrawCaps(caps9.qua, caps9, 60, BLACK);
	Capsule caps10 = { {x - boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },{x - boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z + boxrounded.box.demSize.z },boxrounded.rayon,QuaternionFromAxisAngle({ 1,0,0 }, PI / 2) };//BttomLeft et BottomRight
	MyDrawCaps(caps10.qua, caps10, 60, BLACK);

	//Left face
	Capsule caps11 = { {x + boxrounded.box.demSize.x , y + boxrounded.box.demSize.y , z + boxrounded.box.demSize.z },{x + boxrounded.box.demSize.x , y + boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },boxrounded.rayon ,QuaternionFromAxisAngle({ 1,0,0 }, -PI / 2) };//TopLeft et TopRight
	MyDrawCaps(caps11.qua, caps11, 60, BLACK);
	Capsule caps12 = { {x + boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z + boxrounded.box.demSize.z },{x + boxrounded.box.demSize.x , y - boxrounded.box.demSize.y , z - boxrounded.box.demSize.z },boxrounded.rayon,QuaternionFromAxisAngle({ 1,0,0 }, -PI / 2) };//BottomRight et TopRight
	MyDrawCaps(caps12.qua, caps12, 60, BLACK);

	rlPopMatrix();
}

//Calcul des différentes interesetions dont on a besoin
bool InterSegPlane(Plane plan, Segment segment, Vector3& interPt, Vector3& interNormal) {
	Vector3 AB = Vector3Subtract(segment.pt2, segment.pt1);
	float dotABn = Vector3DotProduct(AB, plan.normal);
	if (fabs(dotABn) < EPSILON) {
		return false;
	}
	float t = (plan.d - Vector3DotProduct(segment.pt1, plan.normal)) / dotABn;
	if (t < 0 || t>1) {
		return false;
	}
	interPt = Vector3Add(segment.pt1, Vector3Scale(AB, t));
	if (dotABn < 0) {
		interNormal = plan.normal;
	}
	else {
		interNormal = Vector3Scale(plan.normal, -1);//Prends -normal
	}
	return true;

}

bool InterSegmentSphere(Segment seg, Sphere s, Vector3& interPt, Vector3& interNormal) {
	Vector3 AB = Vector3Subtract(seg.pt2, seg.pt1);
	Vector3 omegaA = Vector3Subtract(seg.pt1, s.center);
	float delta = 4 * pow(Vector3DotProduct(AB, omegaA), 2) + 4 * Vector3LengthSqr(AB) * pow(s.r, 2) - 4 * Vector3LengthSqr(AB) * Vector3LengthSqr(omegaA);
	float t1 = (-2 * Vector3DotProduct(AB, omegaA) + sqrt(delta)) / (2 * Vector3LengthSqr(AB));
	float t2 = (-2 * Vector3DotProduct(AB, omegaA) - sqrt(delta)) / (2 * Vector3LengthSqr(AB));
	float t = 0;
	if (delta < -EPSILON) {
		return false;
	}
	else if (delta == 0) {
		//t1=t2
		t = t1;
	}
	else {
		if (t2 >= 0 && t2 <= 1) {//t2 est la plus petite valeur des deux
			t = t2;
		}
		else {
			return false;
		}
	}
	interPt = Vector3Add(seg.pt1, Vector3Scale(AB, t));
	Vector3 vecInterPt_centre = Vector3Subtract(interPt, s.center);
	interNormal = Vector3Normalize(vecInterPt_centre);
	return true;
}

bool InterSegDisc(Segment seg, Disc disc, Vector3& interPt, Vector3& interNormal) {
	Plane plan = { disc.ref.j, Vector3DotProduct(disc.ref.j,disc.ref.origin) };
	if (InterSegPlane(plan, seg, interPt, interNormal) == false) {
		return false;
	}
	Vector3 oPrimeI = Vector3Subtract(interPt, disc.ref.origin);
	if (Vector3LengthSqr(oPrimeI) > pow(disc.radius, 2)) {//Car Vector3LengthSqr(oPrimeI) renvoie la norme au carré
		return false;
	}
	return true;
}

bool InterSegmentInfiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal) {
	Vector3 PQ = Vector3Subtract(cyl.pt2, cyl.pt1);
	Vector3 PA = Vector3Subtract(seg.pt1, cyl.pt1);
	Vector3 AB = Vector3Subtract(seg.pt2, seg.pt1);
	Vector3 temp1 = Vector3Subtract(AB, Vector3Scale(PQ, Vector3DotProduct(AB, PQ) / Vector3LengthSqr(PQ)));
	Vector3 temp2 = Vector3Subtract(PA, Vector3Scale(PQ, Vector3DotProduct(PA, PQ) / Vector3LengthSqr(PQ)));

	float a = Vector3LengthSqr(temp1);
	float b = 2 * Vector3DotProduct(temp1, temp2);
	float c = Vector3LengthSqr(temp2) - pow(cyl.radius, 2);

	float delta = b * b - 4 * a * c;
	float t1 = (-b + sqrt(delta)) / (2 * a);
	float t2 = (-b - sqrt(delta)) / (2 * a);
	float t = 0;
	if (delta < 0) {
		printf("delta<0");
		return false;
	}
	else if (delta == 0) {
		//t1=t2
		t = t1;
	}
	else {
		if (t2 >= 0 && t2 <= 1) {//t2 est la plus petite valeur des deux
			t = t2;
		}
		else {
			printf("Il ya un soucis dans le calcul");
		}
	}
	
	interPt = Vector3Add(seg.pt1, Vector3Scale(AB, t));
	//Calcul de la normale
	//On calcule la distance du centre de cylindre au point d'intersection.
	Vector3 vecInterPt_centre = Vector3Subtract(interPt, cyl.pt1);
	Vector3 u = Vector3Normalize(PQ);//On normalise l'axe du cylindre, puis on effectue le projeté orthogonal de p
	float vecInterPt_centre_u = Vector3DotProduct(vecInterPt_centre, u);//On effectue le projeté orthogonal de pi sur l'axe u
	Vector3 ph = Vector3Scale(u, vecInterPt_centre_u);//On a la distance 
	interNormal = Vector3Normalize(Vector3Subtract(interPt, Vector3Add(cyl.pt1, ph)));//on trouve le vecteur normal 
	return true;
}
bool InterSegmentFiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal) {

	Vector3 PQ = Vector3Subtract(cyl.pt2, cyl.pt1);
	Vector3 PA = Vector3Subtract(seg.pt1, cyl.pt1);
	Vector3 AB = Vector3Subtract(seg.pt2, seg.pt1);
	Vector3 u = Vector3Normalize(PQ);
	Vector3 vecInterPt_centre=PA;//on initialise les valeur que l'on recouvrira ensuite. 
	interPt = seg.pt1;//on initialise les valeur que l'on recouvrira ensuite. 


	if (Vector3LengthSqr(Vector3CrossProduct(PA, u)) > pow(cyl.radius, 2)) { //on vérifie que le point est à l'extérieur du cylindre cf:https://fr.wikipedia.org/wiki/Distance_d%27un_point_à_une_droite
		Vector3 temp1 = Vector3Subtract(AB, Vector3Scale(PQ, Vector3DotProduct(AB, PQ) / Vector3LengthSqr(PQ)));
		Vector3 temp2 = Vector3Subtract(PA, Vector3Scale(PQ, Vector3DotProduct(PA, PQ) / Vector3LengthSqr(PQ)));

		float a = Vector3LengthSqr(temp1);
		float b = 2 * Vector3DotProduct(temp1, temp2);
		float c = Vector3LengthSqr(temp2) - pow(cyl.radius, 2);

		float delta = b * b - 4 * a * c;
		float t1 = (-b + sqrt(delta)) / (2 * a);
		float t2 = (-b - sqrt(delta)) / (2 * a);
		float t = 0;
		if (delta < 0) {
			return false;
		}
		else if (delta == 0) {
			//t1=t2
			t = t1;
		}
		else {
			if (t2 >= 0 && t2 <= 1) {//t2 est la plus petite valeur des deux
				t = t2;
			}
			else {
				printf("Il ya un soucis dans le calcul");
			}
		}

		interPt = Vector3Add(seg.pt1, Vector3Scale(AB, t));
		//On calcule la distance du centre de cylindre au point d'intersection.
		Vector3 vecInterPt_centre = Vector3Subtract(interPt, cyl.pt1);
	}

	float vecInterPt_centre_pq = Vector3DotProduct(vecInterPt_centre, PQ);
	Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), 0);
	Vector3 i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
	Vector3 j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	Vector3 k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	if (vecInterPt_centre_pq < 0) { //Au dela de la capsule arrière, en dessous du cylindre, le produit scalaire sera négatif car l'angle sera compris entre PI/2 de PI modulo PI 

		return InterSegDisc(seg, { { {0,2,0},i,j,k },cyl.radius }, interPt, interNormal);
	}
	else if (vecInterPt_centre_pq/Vector3Length(PQ) > Vector3Length(PQ)) { //Au dela de la capsule avant, au dessus du cylindre
		return InterSegDisc(seg, { { {0,2,0},i,j,k },cyl.radius }, interPt, interNormal);
	}
	else { //Dans le cylindre seulement si le commencement est à l'extérieur
		float vecInterPt_centre_u = Vector3DotProduct(vecInterPt_centre, u);
		Vector3 ph = Vector3Scale(u, vecInterPt_centre_u);
		interNormal = Vector3Normalize(Vector3Subtract(interPt, Vector3Add(cyl.pt1, ph)));
		return true;
	}
}

bool InterSegmentQuad(Segment seg, Quad quad, Vector3& interPt, Vector3& interNormal) {
	Plane plan = { quad.ref.j, Vector3DotProduct(quad.ref.j,quad.ref.origin) };
	if (InterSegPlane(plan, seg, interPt, interNormal) == false) {
		return false;
	}
	Vector3 vec = GlobalToLocalPos(interPt, quad.ref);
	if (fabs(vec.x) > quad.demSize.x || fabs(vec.z) > quad.demSize.y) {
		return false;
	}
	return true;
}

bool InterSegmentCapsule(Segment seg, Capsule capsule, Vector3& interPt, Vector3& interNormal) {

	Vector3 PQ = Vector3Subtract(capsule.pt2, capsule.pt1);
	Vector3 PA = Vector3Subtract(seg.pt1, capsule.pt1);
	Vector3 AB = Vector3Subtract(seg.pt2, seg.pt1);
	Vector3 u = Vector3Normalize(PQ);
	Vector3 vecInterPt_centre = { 0 };//on initialise les valeur que l'on recouvrira ensuite. 
	interPt = seg.pt1;//on initialise les valeur que l'on recouvrira ensuite. 


	if (Vector3LengthSqr(Vector3CrossProduct(PA, u)) > pow(capsule.radius, 2)) { //on vérifie que le point est à l'extérieur du cylindre cf:https://fr.wikipedia.org/wiki/Distance_d%27un_point_à_une_droite
		Vector3 temp1 = Vector3Subtract(AB, Vector3Scale(PQ, Vector3DotProduct(AB, PQ) / Vector3LengthSqr(PQ)));
		Vector3 temp2 = Vector3Subtract(PA, Vector3Scale(PQ, Vector3DotProduct(PA, PQ) / Vector3LengthSqr(PQ)));

		float a = Vector3LengthSqr(temp1);
		float b = 2 * Vector3DotProduct(temp1, temp2);
		float c = Vector3LengthSqr(temp2) - pow(capsule.radius, 2);

		float delta = pow(b,2) - 4 * a * c;
		float t1 = (-b + sqrt(delta)) / (2 * a);
		float t2 = (-b - sqrt(delta)) / (2 * a);
		float t = 0;
		if (delta < -EPSILON) {
			//printf("delta<0");
			return false;
		}
		if(delta<EPSILON){
			t = -b / (2 * a);
		}
		else {
			//t2 est la plus petite valeur des deux
			t = t2;
		}
		if (t < -EPSILON || t > 1 + EPSILON) {
			return false;
		}
		interPt = Vector3Add(seg.pt1, Vector3Scale(AB, t));
		//On calcule la distance du centre de cylindre au point d'intersection.
		vecInterPt_centre = Vector3Subtract(interPt, capsule.pt1);
	}
	

	float vecInterPt_centre_pq = Vector3DotProduct(vecInterPt_centre, PQ);
	if (vecInterPt_centre_pq < EPSILON) { //Au dela de la capsule arrière, en dessous du cylindre, le produit scalaire sera négatif car l'angle sera compris entre PI/2 et PI modulo PI 
		return InterSegmentSphere(seg, { capsule.pt1,capsule.radius }, interPt, interNormal);
	}
	else if (vecInterPt_centre_pq > Vector3LengthSqr(PQ) - EPSILON) { //Au dela de la capsule avant, au dessus du cylindre
		return InterSegmentSphere(seg, { capsule.pt2,capsule.radius }, interPt, interNormal);
	}
	else { //Dans le cylindre seulement possible si le commencement est en dehors du cylindre
		float vecInterPt_centre_u = Vector3DotProduct(vecInterPt_centre, u);
		Vector3 ph = Vector3Scale(u, vecInterPt_centre_u);
		interNormal = Vector3Normalize(Vector3Subtract(interPt, Vector3Add(capsule.pt1, ph)));
		return true;
	}
}

bool IntersectSegmentBoxRounded(Segment segment, Box_arrondie box, Vector3& interPt, Vector3& interNormal) {
	Quad * listQuad = ListQuadrilatere(box);
	float dist = -5;
	Vector3 interPtTest;
	Vector3 interNormalTest;
	Vector3 interPtMin;
	Vector3 interNormalMin;
	for (int i = 0; i < 6; i++) {
		if (InterSegmentQuad(segment, listQuad[i], interPtTest, interNormalTest)==true) {
			//On doit choisir le point d'intersection le plus proche pour l'ensemble de TOUS les quads et l'ensemble de TOUTES les Capsules
			//Pareil pour la normale
			float distQuad = Vector3LengthSqr(Vector3Subtract(interPtTest, segment.pt1));
			if (dist<0||distQuad < dist) {
				dist = distQuad;
				interPtMin = interPtTest;
				interNormalMin = interNormalTest;
			}
		}
	}
	Capsule* listCaps = ListCapsule(box);
	for (int i = 0; i < 12; i++) {
		if (InterSegmentCapsule(segment, listCaps[i], interPtTest, interNormalTest)==true) {
			float distCaps = Vector3LengthSqr(Vector3Subtract(interPtTest, segment.pt1));
			if (dist < 0 || (distCaps > 0 && distCaps < dist)) {//Pour etre sur que dist != 0
				dist = distCaps;
				interPtMin = interPtTest;
				interNormalMin = interNormalTest;
			}
		}
	}
	if (dist> 0){//S'il y a eu une intersection
		interPt = interPtMin;
		interNormal = interNormalMin;
		return true;
	}
	return false;
}

bool Collision(Vector3 a, Vector3 b, Box_arrondie * tabBox, Ball &balle, int loop) {
	//On calcule le pt d'arrivé hypothétique, création du segment AB
	//calcul des intersections pour avoir I
	//Calcul du nouveau vecteur vitesse avec reflect
	//Calcul du relicat de distance
	//On positionne la balle
	int lastCollide;
	int nb = 11;
	do {
		Segment seg = { a, b };
		lastCollide = -1;
		for (int i = 0; i < nb; i++) {
			Vector3 interPt;
			Vector3 interNormal;
			if (i != lastCollide && IntersectSegmentBoxRounded(seg, tabBox[i], interPt, interNormal) == true) {
				balle.mvmt = Vector3Reflect(balle.mvmt, interNormal);
				Vector3 relicat = Vector3Subtract(interPt, b);
				balle.position = Vector3Add(interPt, Vector3Scale(Vector3Normalize(balle.mvmt), Vector3Length(relicat)));
				lastCollide = i;
				return true;
			}
		}
		balle.position = b;
		return false;
	} while (lastCollide >= 0);
		
}

int main(int argc, char* argv[])
{
	// Initialization
	//--------------------------------------------------------------------------------------
	float screenSizeCoef = .7f;
	const int screenWidth = 1920 * screenSizeCoef;
	const int screenHeight = 1080 * screenSizeCoef;

	InitWindow(screenWidth, screenHeight, "Gaby's intersection");

	SetTargetFPS(60);

	//CAMERA
	Vector3 cameraPos = { 8.0f, 15.0f, 14.0f };
	Camera camera = { 0 };
	camera.position = cameraPos;
	camera.target = { 0.0f, 0.0f, 0.0f };
	camera.up = { 0.0f, 1.0f, 0.0f };
	camera.fovy = 45.0f;
	camera.type = CAMERA_PERSPECTIVE;
	SetCameraMode(camera, CAMERA_CUSTOM);  // Set an orbital camera mode

	//Sert pour la collsion 
	Ball balle = { 1, { 0,6,0 }, {0,-1,0} };
	Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI / 4);
	Vector3 i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
	Vector3 j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	Vector3 k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	Referential ref = Referential{ {0,-5,0},i,j,k ,qOrient };

	Box_arrondie obstacles[11];
	//Definition de l'arène
	Box box = { ref,{1,1,1} };
	Box_arrondie boxrounded0 = { balle.rayon, box, NULL, NULL };
	obstacles[0] = boxrounded0;
	qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), PI/2);//Gauche
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {0,0,10},i,j,k ,qOrient };
	Box box1= { ref,{10,0,10} };
	Box_arrondie boxrounded1 = { balle.rayon, box1, NULL, NULL };
	obstacles[1] = boxrounded1;
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);//Droite
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {0,0,-10},i,j,k ,qOrient };
	Box box2 = { ref,{10,0,10} };
	Box_arrondie boxrounded2 = { balle.rayon, box2, NULL, NULL };
	obstacles[2] = boxrounded2;
	qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), 0);//Au dessus
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {0,10,0},i,j,k ,qOrient };
	Box box3 = { ref,{10,0,10} };
	Box_arrondie boxrounded3 = { balle.rayon, box3, NULL, NULL };
	obstacles[3] = boxrounded3;
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);//Au dessous
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {0,-10,0},i,j,k ,qOrient };
	Box box4 = { ref,{10,0,10} };
	Box_arrondie boxrounded4 = { balle.rayon, box4, NULL, NULL };
	obstacles[4] = boxrounded4;
	qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 2);//Devant
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {10,0,0},i,j,k ,qOrient };
	Box box5 = { ref,{10,0,10} };
	Box_arrondie boxrounded5 = { balle.rayon, box5, NULL, NULL };
	obstacles[5] = boxrounded5;
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);//Derriere
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {-10,0,0},i,j,k ,qOrient };
	Box box6 = { ref,{10,0,10} };
	Box_arrondie boxrounded6 = { balle.rayon, box6, NULL, NULL };
	obstacles[6] = boxrounded6;
	//Boxs test
	qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI);
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {4,4,4},i,j,k ,qOrient };
	Box box7 = { ref,{3.2f,1,2} };
	Box_arrondie boxrounded7 = { balle.rayon, box7, NULL, NULL };
	obstacles[7] = boxrounded7;
	qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 8);
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {6,3,-2},i,j,k ,qOrient };
	Box box8 = { ref,{1,2,1.5f} };
	Box_arrondie boxrounded8 = { balle.rayon, box8, NULL, NULL };
	obstacles[8] = boxrounded8;
	qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 3);
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {-4,1,1},i,j,k ,qOrient };
	Box box9 = { ref,{4,1,1} };
	Box_arrondie boxrounded9 = { balle.rayon, box9, NULL, NULL };
	obstacles[9] = boxrounded9;
	qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 0,0,1 }), PI / 5);//Devant
	i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
	j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
	k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
	ref = Referential{ {6,-5,8},i,j,k ,qOrient };
	ref = Referential{ {6,-5,8},i,j,k ,qOrient };
	Box box10 = { ref,{1,5,1.8f} };
	Box_arrondie boxrounded10 = { balle.rayon, box10, NULL, NULL };
	obstacles[10] = boxrounded10;

	//--------------------------------------------------------------------------------------

	// Main game loop
	while (!WindowShouldClose())    // Detect window close button or ESC key
	{
		// Update
		//----------------------------------------------------------------------------------
		// TODO: Update your variables here
		//----------------------------------------------------------------------------------

		float deltaTime = GetFrameTime();
		float time = (float)GetTime();

		MyUpdateOrbitalCamera(&camera, deltaTime);

		// Draw
		//----------------------------------------------------------------------------------
		BeginDrawing();

		ClearBackground(RAYWHITE);

		BeginMode3D(camera);
		{
			//Test implémentation quaternion 
			/*Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({1,3,-4}), time);
			MyDrawSphereEx2(qOrient, Vector3{ 0 }, Vector3{2f,2f,2f}, 40, 20, BLUE);
			MyDrawSphereWiresEx2(qOrient, Vector3{ 0 }, Vector3{2f,2f,2f}, 40, 20, WHITE);*/

			//Test fonction demandées intersection segment plan;
			/*Plane p = {Vector3{0.0f,1.0f,0.0f},1.5f};
			Segment s = { Vector3{0.0f,4.0f,0.0f},Vector3{3.0f,0.0f,0.0f} };
			Vector3 interPt;
			Vector3 interNormal;
			InterSegPlane(p, s, interPt, interNormal);
			MyDrawQuad(Vector3Scale(p.normal, p.d), Vector2{ 5.0f,5.0f }, GRAY);
			MyDrawQuadWire(Vector3Scale(p.normal, p.d), Vector2{ 5.0f,5.0f }, DARKGRAY);
			DrawLine3D(s.pt1, s.pt2, RED);
			Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), time);
			MyDrawSphereEx2(qOrient, interPt, Vector3{0.2f,0.2f,0.2f}, 40, 20, BLUE);*/

			//Test Draw Cylinder
			//Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), 45);
			//MyDrawDisk(qOrient, Vector3{3,3,3}, 2, 40, GREEN);
			//MyDrawCylinder(qOrient, {{1,1,1},{2,3,0},1}, 60, true, YELLOW);

			//Test Intersection segment / sphere
			/*Segment seg = {Vector3{0.0f,4.0f,0.0f},Vector3{3.0f,0.0f,0.0f}};
			Sphere s = { Vector3{2.0f,2.0f,0.0f}, Vector3{1.2f,1.2f,1.2f}};
			Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 0,0,0 }), time);
			Vector3 interPt;
			Vector3 interNormal;
			InterSegmentSphere(seg, s, interPt, interNormal);
			DrawLine3D(seg.pt1, seg.pt2, RED);
			MyDrawSphereWiresEx2(qOrient, s.center, s.r, 40, 20, Color{ 100, 0, 50, 180 });
			MyDrawSphereEx2(qOrient, seg.pt1, Vector3{0.2f,0.2f,0.2f}, 40, 20, RED);
			MyDrawSphereEx2(qOrient, seg.pt2, Vector3{0.2f,0.2f,0.2f}, 40, 20, GREEN);
			MyDrawSphereEx2(qOrient, interPt, Vector3{0.2f,0.2f,0.2f}, 40, 20, BLUE);*/

			//Test Intersection segment / infCylinder -> MARCHE
			/*Segment seg = {Vector3{0.0f,4.0f,0.0f},Vector3{3.0f,0.0f,0.0f}};
			Cylinder cyl = { {2,0,0},{2,3,0},0.5f };
			Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), 0);
			Vector3 interPt;
			Vector3 interNormal;
			InterSegmentInfiniteCylinder(seg, cyl, interPt,interNormal);
			DrawLine3D(seg.pt1, seg.pt2, RED);
			MyDrawCylinder(qOrient, cyl, 40, true, YELLOW);
			MyDrawSphereEx2(qOrient, seg.pt1, Vector3{0.2f,0.2f,0.2f}, 40, 20, RED);
			MyDrawSphereEx2(qOrient, seg.pt2, Vector3{0.2f,0.2f,0.2f}, 40, 20, GREEN);
			MyDrawSphereEx2(qOrient, interPt, Vector3{0.2f,0.2f,0.2f}, 40, 20, BLUE);*/

			//Test Intersection segment / Cylinder finis-> MARCHE
			/*Segment seg = {Vector3{0.0f,4.0f,0.0f},Vector3{3.0f,0.0f,0.0f}};
			Cylinder cyl = { {2,0,0},{2,3,0},0.5f };
			Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), 0);
			Vector3 interPt;
			Vector3 interNormal;
			InterSegmentFiniteCylinder(seg, cyl, interPt, interNormal);
			DrawLine3D(seg.pt1, seg.pt2, RED);
			MyDrawCylinder(qOrient, cyl, 40, true, YELLOW);
			MyDrawSphereEx2(qOrient, seg.pt1, Vector3{ 0.2f,0.2f,0.2f }, 40, 20, RED);
			MyDrawSphereEx2(qOrient, seg.pt2, Vector3{ 0.2f,0.2f,0.2f }, 40, 20, GREEN);
			MyDrawSphereEx2(qOrient, interPt, Vector3{ 0.2f,0.2f,0.2f }, 40, 20, BLUE);*/

			//Test Intersection segment / Capsule -> MARCHE
			/*Segment seg = {Vector3{0.0f,4.0f,0.0f},Vector3{3.0f,0.0f,0.0f}};
			Capsule caps = { {2,0,0},{2,3,0},0.5 };
			Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), 0);
			Vector3 interPt;
			Vector3 interNormal;
			InterSegmentCapsule(seg, caps, interPt, interNormal);
			DrawLine3D(seg.pt1, seg.pt2, RED);
			MyDrawCaps(qOrient, caps, 40, GREEN);
			MyDrawSphereEx2(qOrient, seg.pt1, Vector3{ 0.2f,0.2f,0.2f }, 40, 20, RED);
			MyDrawSphereEx2(qOrient, seg.pt2, Vector3{ 0.2f,0.2f,0.2f }, 40, 20, GREEN);
			MyDrawSphereEx2(qOrient, interPt, Vector3{ 0.2f,0.2f,0.2f }, 40, 20, BLUE);*/

			//Test InterSegDisc -> MARCHE
			/*Segment seg = {Vector3{0.0f,4.0f,0.0f},Vector3{3.0f,0.0f,0.0f}};
			Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), time);
			Vector3 i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
			Vector3 j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
			Vector3 k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
			Referential ref = Referential{ {0,2,0},i,j,k };
			Disc disk = { ref,3};
			Vector3 interPt;
			Vector3 interNormal;
			InterSegDisc(seg, disk, interPt, interNormal);
			DrawLine3D(seg.pt1, seg.pt2, RED);
			MyDrawDisk(qOrient, disk.ref.origin, disk.radius, 40, YELLOW);
			MyDrawSphereEx2(qOrient, seg.pt1, Vector3{0.2f,0.2f,0.2f}, 40, 20, RED);
			MyDrawSphereEx2(qOrient, seg.pt2, Vector3{0.2f,0.2f,0.2f}, 40, 20, GREEN);
			MyDrawSphereEx2(qOrient, interPt, Vector3{0.2f,0.2f,0.2f} ,40, 20, BLUE);*/

			//Test inter SegQuad -> MARCHE
			/*Segment seg = {Vector3{0.0f,4.0f,0.0f},Vector3{3.0f,0.0f,0.0f}};
			Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), 45);
			Vector3 i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
			Vector3 j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
			Vector3 k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
			Referential ref = Referential{ {0,2,0},i,j,k };
			Quad quad = { ref,{4,0,4} };
			Vector3 interPt;
			Vector3 interNormal;
			InterSegmentQuad(seg, quad, interPt, interNormal);
			DrawLine3D(seg.pt1, seg.pt2, RED);
			MyDrawQuad(qOrient, quad.ref.origin, { quad.demSize.x, quad.demSize.z }, YELLOW);
			MyDrawSphereEx2(qOrient, seg.pt1, Vector3{0.2f,0.2f,0.2f}, 40, 20, RED);
			MyDrawSphereEx2(qOrient, seg.pt2, Vector3{0.2f,0.2f,0.2f}, 40, 20, GREEN);
			MyDrawSphereEx2(qOrient, interPt, Vector3{0.2f,0.2f,0.2f}, 40, 20, BLUE);*/

			//Test Draw Caps
			/*Capsule caps = {{0,0,0},{0,10,0},0.5};
			Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({ 1,0,0 }), 0);
			MyDrawCaps(qOrient, caps, 40, GREEN);*/

			//Test MyDrawBox
			/*Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({1,0,0}), 0);
			Vector3 i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
			Vector3 j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
			Vector3 k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
			Referential ref = Referential{ {1,1,1},i,j,k, qOrient };
			Box box = { ref,{2,1,4} };
			MyDrawBox(qOrient,box, BLUE);*/

			//Test MyDrawBoundingBox
			/*Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({1,0,0}), PI / 4);
			Vector3 i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
			Vector3 j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
			Vector3 k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
			Referential ref = Referential{ {1,0,0},i,j,k,qOrient};
			Box box = { ref,{2,1,4} };
			Box_arrondie boxrounded = {0.3f, box, NULL, NULL};
			MyDrawRoundedBox(qOrient,boxrounded, ORANGE);*/

			//test du listing d'élément de la rounded box
			/*Quad* tab = ListQuadrilatere(boxrounded);
			MyDrawQuad(tab[0].ref.qua, tab[0].ref.origin, tab[0].demSize, RED);
			MyDrawQuad(tab[1].ref.qua, tab[1].ref.origin, tab[1].demSize, RED);
			MyDrawQuad(tab[2].ref.qua, tab[2].ref.origin, tab[2].demSize, RED);
			MyDrawQuad(tab[3].ref.qua, tab[3].ref.origin, tab[3].demSize, RED);
			MyDrawQuad(tab[4].ref.qua, tab[4].ref.origin, tab[4].demSize, RED);
			MyDrawQuad(tab[5].ref.qua, tab[5].ref.origin, tab[5].demSize, RED);

			Capsule* tabC = ListCapsule(boxrounded);
			MyDrawCaps(tabC[0].qua, tabC[0], 60, PURPLE);
			MyDrawCaps(tabC[1].qua, tabC[1], 60, PURPLE);
			MyDrawCaps(tabC[2].qua, tabC[2], 60, PURPLE);
			MyDrawCaps(tabC[3].qua, tabC[3], 60, PURPLE);
			MyDrawCaps(tabC[4].qua, tabC[4], 60, PURPLE);
			MyDrawCaps(tabC[5].qua, tabC[5], 60, PURPLE);
			MyDrawCaps(tabC[6].qua, tabC[6], 60, PURPLE);
			MyDrawCaps(tabC[7].qua, tabC[7], 60, PURPLE);
			MyDrawCaps(tabC[8].qua, tabC[8], 60, PURPLE);
			MyDrawCaps(tabC[9].qua, tabC[9], 60, PURPLE);
			MyDrawCaps(tabC[10].qua, tabC[10], 60, PURPLE);
			MyDrawCaps(tabC[11].qua, tabC[11], 60, PURPLE);*/
		
			//Test Intersection segment / box rounded -> MARCHE 
			/*Segment seg = {Vector3{0.0f,4.0f,0.0f},Vector3{4.5f,0.0f,0.0f}};
			Quaternion qOrient = QuaternionFromAxisAngle(Vector3Normalize({1,0,0}), 0);
			Vector3 i = Vector3RotateByQuaternion({ 1, 0, 0 }, qOrient);
			Vector3 j = Vector3RotateByQuaternion({ 0, 1, 0 }, qOrient);
			Vector3 k = Vector3RotateByQuaternion({ 0, 0, 1 }, qOrient);
			Referential ref = Referential{ {1,0,0},i,j,k ,qOrient};
			Box box = { ref,{2,1,4} };
			Box_arrondie boxrounded = {0.3f, box, NULL, NULL};
			Vector3 interPt;
			Vector3 interNormal;
			IntersectSegmentBoxRounded(seg, boxrounded, interPt, interNormal);
			DrawLine3D(seg.pt1, seg.pt2, RED);
			MyDrawRoundedBox(qOrient, boxrounded, GREEN);
			MyDrawSphereEx2(qOrient, seg.pt1, Vector3{ 0.2f,0.2f,0.2f }, 40, 20, RED);
			MyDrawSphereEx2(qOrient, seg.pt2, Vector3{ 0.2f,0.2f,0.2f }, 40, 20, GREEN);
			MyDrawSphereEx2(qOrient, interPt, Vector3{ 0.2f,0.2f,0.2f }, 40, 20, BLUE);*/

			
			if (deltaTime > 0) {
				MyDrawSphereEx2(qOrient, balle.position, { balle.rayon,balle.rayon,balle.rayon }, 60, 60, DARKPURPLE);
				MyDrawSphereWiresEx2(qOrient, balle.position, { balle.rayon,balle.rayon,balle.rayon }, 60, 60, WHITE);
				//Box début
				MyDrawBox(qOrient, box, DARKBLUE);
				//ARENE
				MyDrawBox(qOrient, box1, ORANGE);
				MyDrawWiresBox(qOrient, box1, BLACK);
				MyDrawBox(qOrient, box2, ORANGE);
				MyDrawWiresBox(qOrient, box2, BLACK);
				//MyDrawBox(qOrient, box3, ORANGE);
				MyDrawBox(qOrient, box4, ORANGE);
				MyDrawWiresBox(qOrient, box4, BLACK);
				MyDrawBox(qOrient, box5, ORANGE);
				MyDrawWiresBox(qOrient, box5, BLACK);
				//MyDrawBox(qOrient, box6, ORANGE);
				//Test total
				MyDrawBox(qOrient, box7, MAGENTA);
				MyDrawBox(qOrient, box8, RED);
				MyDrawBox(qOrient, box9, GREEN);
				MyDrawBox(qOrient, box10, DARKGREEN);
				balle.mvmt = Vector3Add(balle.mvmt, Vector3Scale(GRAVITY, deltaTime));//On mets la gravity
				Vector3 b = Vector3Add(balle.position, Vector3Scale(balle.mvmt, deltaTime));//Calcul du point d'arrivée
				bool collide = Collision(balle.position, b, obstacles, balle, deltaTime);
				//printf("%d", collide);
			}
			

			//3D REFERENTIAL
			/*DrawGrid(20, 1.0f);        // Draw a grid
			DrawLine3D({ 0 }, { 0,10,0 }, DARKGRAY);
			DrawSphere({ 10,0,0 }, .2f, RED);
			DrawSphere({ 0,10,0 }, .2f, GREEN);
			DrawSphere({ 0,0,10 }, .2f, BLUE);*/
		}
		EndMode3D();

		EndDrawing();
		//----------------------------------------------------------------------------------
	}


	// De-Initialization
	//--------------------------------------------------------------------------------------   
	CloseWindow();        // Close window and OpenGL context
	//--------------------------------------------------------------------------------------

	return 0;
}

