#include "PBD.hpp"
using namespace mdx;


bool  PBD::init_ShapeMatchingConstraint(const std::vector<Point3d> x0, const std::vector<double> invMasses,
                                  int numPoints,
                                  Point3d& restCm,
                                  Matrix3d& invRestMat) {
    invRestMat.identity();

    // center of mass
    restCm = 0;
    double wsum = 0.0;
    for(int i = 0; i < numPoints; i++) {
        double wi = 1. / (invMasses[i] + EPS);
        restCm += x0[i] * wi;
        wsum += wi;
    }
    if(wsum == 0.0)
        return false;
    restCm /= wsum;

    // A
    Matrix3d A;
    A = 0;
    for(int i = 0; i < numPoints; i++) {
        const Point3d qi = x0[i] - restCm;
        double wi = 1.0 / (invMasses[i] + EPS);
        double x2 = wi * qi[0] * qi[0];
        double y2 = wi * qi[1] * qi[1];
        //double z2 = wi * qi[2] * qi[2];
        double xy = wi * qi[0] * qi[1];
        double xz = wi * qi[0] * qi[2];
        double yz = wi * qi[1] * qi[2];
        A[0] += Point3d(x2, xy, xz);
        A[1] += Point3d(xy, y2, yz);
        //A[2] += Point3d(xz, yz, z2);
        A[2] = Point3d(0, 0, 1);
        ;
    }
    double determinant = det(A);
    if(fabs(determinant) > EPS) {
        invRestMat = inverse(A);
        return true;
    }
    return false;
}

bool  PBD::solve_ShapeMatchingConstraint(const std::vector<Point3d> x0,
                                   const std::vector<Point3d> x,
                                   const std::vector<double> invMasses,
                                   int numPoints,
                                   const Point3d& restCm,
                                   const Matrix3d& invRestMat,
                                   bool allowStretch,
                                   std::vector<Point3d>& corr,
                                   Matrix3d* rot=NULL) {
    if(numPoints <= 0) {
        mdxInfo << "PBD::solve_ShapeMatchingConstraint: numPoints = " << numPoints << endl;
        return false;
    }
    if(restCm == Point3d(0,0,0)) {
        mdxInfo << "PBD::solve_ShapeMatchingConstraint: restCm = " << restCm << endl;
        return false;
    }

    // center of mass
    Point3d cm(0.0, 0.0, 0.0);
    double wsum = 0.0;
    for(int i = 0; i < numPoints; i++) {
        double wi = 1.0 / (invMasses[i] + EPS);
        cm += x[i] * wi;
        wsum += wi;
    }
    if(wsum == 0.0)
        return false;
    cm /= wsum;

    // A
    Matrix3d mat;
    for(int i = 0; i < numPoints; i++) {
        Point3d q = x0[i] - restCm;
        Point3d p = x[i] - cm;
        double w = 1.0 / (invMasses[i] + EPS);
        p *= w;
        mat[0] += Point3d(p[0] * q[0], p[0] * q[1], p[0] * q[2]);
        mat[1] += Point3d(p[1] * q[0], p[1] * q[1], p[1] * q[2]);
        mat[2] += Point3d(p[2] * q[0], p[2] * q[1], p[2] * q[2]);
    }

    mat = mat * invRestMat;
    Matrix3d R, U;
    if (allowStretch)
      R = mat;
    else
      PolarDecompX(mat, R, U);

    for(int i = 0; i < numPoints; i++) {
        Point3d goal = cm + R * (x0[i] - restCm);
        corr.push_back((goal - x[i]));
    }
    if(rot)
        *rot = R;

    return true;
}


bool  PBD::solve_DistanceConstraint(
    const Point3d &p0, double invMass0,
    const Point3d &p1, double invMass1,
    const double restLength,
    const double compressionStiffness,
    const double stretchStiffness,
    Point3d &corr0, Point3d &corr1)
{

    double wSum = invMass0 + invMass1;
    if (wSum == 0.0)
        return false;
    Point3d n = p1 - p0;
    double d = n.norm();
    n /= n.norm();
    Point3d corr;
    if (d < restLength)
        corr = compressionStiffness * n * (d - restLength) / wSum;
    else
        corr = stretchStiffness * n * (d - restLength) / wSum;
    corr0 =  invMass0 * corr;
    corr1 = -invMass1 * corr;
    return true;
}


// XPBD
bool PBD::solve_DistanceConstraint_XPBD(const Point3d &p0, double invMass0,
                               const Point3d &p1, double invMass1,
                               const double restLength,
                               const double compressionStiffness,
                               const double stretchStiffness,
                               Point3d &corr0, Point3d &corr1, double &lambda, double compliance)
        {

    double wSum = invMass0 + invMass1;
    if (wSum == 0.0)
        return false;

    Point3d n = p1 - p0;
    double d = n.norm();
    n /= d;

    double stiffness = 0;
    if (d < restLength)
        stiffness = compressionStiffness;
    else
        stiffness = stretchStiffness;

    Point3d grad = stiffness * n;
    double deltalambda = (-1. * (d - restLength) - (compliance * lambda)) / (wSum + compliance) ;

    lambda += deltalambda;

    corr0 = -deltalambda * invMass0 * grad;
    corr1 = deltalambda * invMass1 * grad;
    return true;

}

// ----------------------------------------------------------------------------------------------
bool  PBD::init_StrainTriangleConstraint(
    const Point3d &p0,
    const Point3d &p1,
    const Point3d &p2,
    Matrix2d &invRestMat)
{
    double a = p1[0] - p0[0]; double b = p2[0] - p0[0];
    double c = p1[1] - p0[1]; double d = p2[1] - p0[1];

    // inverse
    double det = a*d - b*c;
    if (fabs(det) < EPS)
        return false;

    double s = 1. / det;
    invRestMat[0][0] =  d*s;  invRestMat[0][1] = -b*s;
    invRestMat[1][0] = -c*s;  invRestMat[1][1] =  a*s;

    return true;
}

// ----------------------------------------------------------------------------------------------
bool  PBD::solve_StrainTriangleConstraint(
        const Point3d &p0,
        const Point3d &p1,
        const Point3d &p2,
        const Matrix2d &invRestMat,
        const double xxStiffness,
        const double yyStiffness,
        const double xyStiffness,
        const bool normalizeStretch,
        const bool normalizeShear,
        Point3d &corr0, Point3d &corr1, Point3d &corr2, Point3d &lambda, double alpha_tilde)
{
    Point3d c[2];
    c[0] = Point3d(invRestMat[0][0], invRestMat[1][0], 0.0);
    c[1] = Point3d(invRestMat[0][1], invRestMat[1][1], 0.0);

    Point3d r[3];
    //int lambda_i = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j <= i; j++) {

// 			r[0] = Point3d(p1[0] - p0[0], p2[0] - p0[0], 0.0);  // Jacobi
// 			r[1] = Point3d(p1[1] - p0[1], p2[1] - p0[1], 0.0);
// 			r[2] = Point3d(p1[2] - p0[2], p2[2] - p0[2], 0.0);

            r[0] = Point3d((p1[0] + corr1[0]) - (p0[0] + corr0[0]), (p2[0] + corr2[0]) - (p0[0] + corr0[0]), 0.0);		// Gauss - Seidel
            r[1] = Point3d((p1[1] + corr1[1]) - (p0[1] + corr0[1]), (p2[1] + corr2[1]) - (p0[1] + corr0[1]), 0.0);
            r[2] = Point3d((p1[2] + corr1[2]) - (p0[2] + corr0[2]), (p2[2] + corr2[2]) - (p0[2] + corr0[2]), 0.0);


            double Sij = 0.0;
            for (int k = 0; k < 3; k++)
                Sij += (r[k] * c[i]) * (r[k] * c[j]);

            Point3d d[3];
            d[0] = Point3d(0.0, 0.0, 0.0);

            for (int k = 0; k < 2; k++) {
                d[k+1]  = Point3d(r[0] * c[j], r[1] * c[j], r[2] * c[j]) * invRestMat[k][i];
                d[k+1] += Point3d(r[0] * c[i], r[1] * c[i], r[2] * c[i]) * invRestMat[k][j];
                d[0] -= d[k+1];
            }

            if (i != j && normalizeShear) {
                double fi2 = 0.0;
                double fj2 = 0.0;
                for (int k = 0; k < 3; k++) {
                    fi2 += (r[k] * c[i]) * (r[k] * c[i]);
                    fj2 += (r[k] * c[j]) * (r[k] * c[j]);
                }
                double fi = sqrt(fi2);
                double fj = sqrt(fj2);

                d[0] = Point3d(0.0, 0.0, 0.0);
                double s = Sij / (fi2*fi*fj2*fj);
                for (int k = 0; k < 2; k++) {
                    d[k+1] /= fi * fj;
                    d[k+1] -= fj*fj * Point3d(r[0] * (c[i]), r[1] * (c[i]), r[2] * (c[i])) * invRestMat[k][i] * s;
                    d[k+1] -= fi*fi * Point3d(r[0] * (c[j]), r[1] * (c[j]), r[2] * (c[j])) * invRestMat[k][j] * s;
                    d[0] -= d[k+1];
                }
                Sij = Sij / (fi * fj);
            }


            /*
            double numerator = -lambda[lambda_i] * alpha_tilde;
            double denominator =
                 d[0].norm()*d[0].norm() +
                 d[1].norm()*d[1].norm() +
                 d[2].norm()*d[2].norm() + alpha_tilde;
            */

            double lambda =
                 d[0].norm()*d[0].norm() +
                 d[1].norm()*d[1].norm() +
                 d[2].norm()*d[2].norm();
            /*
            if (denominator == 0.0)
                continue;
            */
            if (i == 0 && j == 0) {
                if (normalizeStretch) {
                    double s = sqrt(Sij);
                    lambda = 2. * s * (s - 1.) / lambda * xxStiffness;
                    //numerator += -2. * s * (s - 1.) * xxStiffness;
                }
                else {
                    lambda = (Sij - 1.) / lambda * xxStiffness;
                    //numerator += -(Sij - 1.) * xxStiffness;
                }
            }
            else if (i == 1 && j == 1) {
                if (normalizeStretch) {
                    double s = sqrt(Sij);
                    lambda = 2. * s * (s - 1.) / lambda * yyStiffness;
                    //numerator += -2. * s * (s - 1.) * yyStiffness;
                }
                else {
                    lambda = (Sij - 1.) / lambda * yyStiffness;
                    //numerator += -(Sij - 1.) * yyStiffness;
                }
            }
            else {
                lambda = Sij / lambda * xyStiffness;
                //numerator += -Sij  * xyStiffness;
            }
            //double deltalambda = numerator / denominator;
            corr0 -= lambda  * d[0];
            corr1 -= lambda  * d[1];
            corr2 -= lambda  * d[2];
            //corr0 += deltalambda  * d[0];
            //corr1 += deltalambda  * d[1];
            //corr2 += deltalambda  * d[2];
            //lambda[lambda_i] += deltalambda;
            //lambda_i ++;
        }
    }
    return true;
}

bool  PBD::solve_StrainTriangleConstraint(
         const Point3d &p0, double invMass0,
         const Point3d &p1, double invMass1,
         const Point3d &p2, double invMass2,
         const Matrix2d &invRestMat,
         const double xxStiffness,
         const double yyStiffness,
         const double xyStiffness,
         const bool normalizeStretch,
         const bool normalizeShear,
         Point3d &corr0, Point3d &corr1, Point3d &corr2)
 {
     Point3d c[2];
     c[0] = Point3d(invRestMat[0][0], invRestMat[1][0], 0.0);
     c[1] = Point3d(invRestMat[0][1], invRestMat[1][1], 0.0);

     Point3d r[3];

     for (int i = 0; i < 2; i++) {
         for (int j = 0; j <= i; j++) {

 //          r[0] = Point3d(p1[0] - p0[0], p2[0] - p0[0], 0.0);  // Jacobi
 //          r[1] = Point3d(p1[1] - p0[1], p2[1] - p0[1], 0.0);
 //          r[2] = Point3d(p1[2] - p0[2], p2[2] - p0[2], 0.0);

             r[0] = Point3d((p1[0] + corr1[0]) - (p0[0] + corr0[0]), (p2[0] + corr2[0]) - (p0[0] + corr0[0]), 0.0);     // Gauss - Seidel
             r[1] = Point3d((p1[1] + corr1[1]) - (p0[1] + corr0[1]), (p2[1] + corr2[1]) - (p0[1] + corr0[1]), 0.0);
             r[2] = Point3d((p1[2] + corr1[2]) - (p0[2] + corr0[2]), (p2[2] + corr2[2]) - (p0[2] + corr0[2]), 0.0);


             double Sij = 0.0;
             for (int k = 0; k < 3; k++)
                 Sij += (r[k] * (c[i])) * (r[k] * (c[j]));

             Point3d d[3];
             d[0] = Point3d(0.0, 0.0, 0.0);

             for (int k = 0; k < 2; k++) {
                 d[k+1]  = Point3d((r[0] * (c[j])), (r[1] * (c[j])), (r[2] * (c[j]))) * invRestMat(k, i);
                 d[k+1] += Point3d((r[0] * (c[i])), (r[1] * (c[i])), (r[2] * (c[i]))) * invRestMat(k, j);
                 d[0] -= d[k+1];
             }

             if (i != j && normalizeShear) {
                 double fi2 = 0.0;
                 double fj2 = 0.0;
                 for (int k = 0; k < 3; k++) {
                     fi2 += (r[k] * (c[i])) * (r[k] * (c[i]));
                     fj2 += (r[k] * (c[j])) * (r[k] * (c[j]));
                 }
                 double fi = sqrt(fi2);
                 double fj = sqrt(fj2);

                 d[0] = Point3d(0.0, 0.0, 0.0);
                 double s = Sij / (fi2*fi*fj2*fj);
                 for (int k = 0; k < 2; k++) {
                     d[k+1] /= fi * fj;
                     d[k+1] -= fj*fj * Point3d((r[0] * (c[i])), (r[1] * (c[i])), (r[2] * (c[i]))) * invRestMat[k][i] * s;
                     d[k+1] -= fi*fi * Point3d((r[0] * (c[j])), (r[1] * (c[j])), (r[2] * (c[j]))) * invRestMat[k][j] * s;
                     d[0] -= d[k+1];
                 }
                 Sij = Sij / (fi * fj);
             }

             double lambda =
                 invMass0 * d[0].norm() * d[0].norm()  +
                 invMass1 * d[1].norm() * d[1].norm()  +
                 invMass2 * d[2].norm() * d[2].norm();

             if (lambda == 0.0)
                 continue;

             if (i == 0 && j == 0) {
                 if (normalizeStretch) {
                     double s = sqrt(Sij);
                     lambda = 2. * s * (s - 1.) / lambda * xxStiffness;
                 }
                 else {
                     lambda = (Sij - 1.) / lambda * xxStiffness;
                 }
             }
             else if (i == 1 && j == 1) {
                 if (normalizeStretch) {
                     double s = sqrt(Sij);
                     lambda = 2. * s * (s - 1.) / lambda * yyStiffness;
                 }
                 else {
                     lambda = (Sij - 1.) / lambda * yyStiffness;
                 }
             }
             else {
                 lambda = Sij / lambda * xyStiffness;
             }

             corr0 -= lambda * invMass0 * d[0];
             corr1 -= lambda * invMass1 * d[1];
             corr2 -= lambda * invMass2 * d[2];
         }
     }
     return true;
 }

bool  PBD::solve_PressureConstraint(
       const CCStructure& cs,
       const CCIndexDataAttr &indexAttr,
       Tissue::VertexDataAttr &vMAttr,
       const std::set<CCIndex> fs,
       const double restArea,
       const double pressure,
       std::map<CCIndex, Point3d>& corr, double &lambda, double compliance=0)
{
   std::map<CCIndex, Point3d> grads;
   double area = 0;
   double denominator = 0;
   Point3d p1, p2;
   CCIndex v0, v1, v2;
   for(CCIndex f : fs) {
       area += indexAttr[f].measure;
       std::vector<CCIndex> vsn = faceVertices(cs, f);
       v0 = vsn[0], v1 = vsn[1], v2 = vsn[2];
       p1 = indexAttr[v0].pos - indexAttr[v2].pos;
       p2 = indexAttr[v1].pos - indexAttr[v2].pos;
       grads[v2] += (-2*(p2 ^ (p1 ^ p2)) - 2*(p1 ^ (p2 ^ p1))) /  indexAttr[f].measure;
       p1 = indexAttr[v1].pos - indexAttr[v0].pos;
       p2 = indexAttr[v2].pos - indexAttr[v0].pos;
       grads[v0] += (-2*(p2 ^ (p1 ^ p2)) - 2*(p1 ^ (p2 ^ p1))) /  indexAttr[f].measure;
       p1 = indexAttr[v2].pos - indexAttr[v1].pos;
       p2 = indexAttr[v0].pos - indexAttr[v1].pos;
       grads[v1] += (-2*(p2 ^ (p1 ^ p2)) - 2*(p1 ^ (p2 ^ p1))) /  indexAttr[f].measure;
   }

   for(auto g : grads)
       denominator += vMAttr[g.first].invmass * norm(g.second) * norm(g.second);
   denominator += compliance;
   if(denominator < EPS) return false;
   double constrain = area - restArea*pressure;
   double deltalambda = (-1. * constrain - compliance * lambda) / denominator;
   //double lambda = constrain / denominator;
   lambda += deltalambda;
   for(auto g : grads)
        corr[g.first] = deltalambda * vMAttr[g.first].invmass * g.second;

   return true;
}


bool  PBD::solve_DihedralConstraint(
     const Point3d &p0, double invMass0,
     const Point3d &p1, double invMass1,
     const Point3d &p2, double invMass2,
     const Point3d &p3, double invMass3,
     const double restAngle,
     const double stiffness,
     Point3d &corr0, Point3d &corr1, Point3d &corr2, Point3d &corr3, bool verbose=false)
 {
     // derivatives from Bridson, Simulation of Clothing with Folds and Wrinkles
     // his modes correspond to the derivatives of the bending angle arccos(n1 dot n2) with correct scaling

     if (invMass0 == 0.0 && invMass1 == 0.0)
         return false;

     Point3d e = p3-p2;
     double  elen = e.norm();
     if (elen < EPS)
         return false;

     double invElen = 1.0 / elen;

     Point3d n1 = (p2-p0).cross(p3-p0); n1 /= n1.norm() * n1.norm();
     Point3d n2 = (p3 - p1).cross(p2 - p1); n2 /= n2.norm() * n2.norm();

     Point3d d0 = elen*n1;
     Point3d d1 = elen*n2;
     Point3d d2 = ((p0-p3) * e) * invElen * n1 + ((p1-p3) * e) * invElen * n2;
     Point3d d3 = ((p2-p0) * e) * invElen * n1 + ((p2-p1) * e) * invElen * n2;

     n1 /= n1.norm();
     n2 /= n2.norm();
     double dot = n1 * n2;

     if (dot < -1.0) dot = -1.0;
     if (dot >  1.0) dot =  1.0;
     double phi = acos(dot);

     if(verbose)
         mdxInfo << "Dihedral angle: " << phi* (180./M_PI) << " restAngle: " << restAngle* (180./M_PI)  <<  endl;

     //if(phi<restAngle) return false;

     // double phi = (-0.6981317 * dot * dot - 0.8726646) * dot + 1.570796;    // fast approximation

     double lambda =
         invMass0 * d0.norm() * d0.norm() +
         invMass1 * d1.norm() * d1.norm() +
         invMass2 * d2.norm() * d2.norm() +
         invMass3 * d3.norm() * d3.norm();

     if (lambda == 0.0)
         return false;

     // stability
     // 1.5 is the largest magic number I found to be stable in all cases :-)
     //if (stiffness > 0.5 && fabs(phi - b.restAngle) > 1.5)
     //  stiffness = 0.5;

     lambda = (phi - restAngle) / lambda * stiffness;

     if (n1.cross(n2) * e > 0.0) ///////////// I changed < to >
         lambda = -lambda;

    corr0 = - invMass0 * lambda * d0;
    corr1 = - invMass1 * lambda * d1;
    corr2 = - invMass2 * lambda * d2;
    corr3 = - invMass3 * lambda * d3;

    return true;
}



// ----------------------------------------------------------------------------------------------
bool PBD::init_FEMTriangleConstraint(
        const Point3d &p0,
        const Point3d &p1,
        const Point3d &p2,
        double &area,
        Matrix2d &invRestMat)
{
        Point3d normal0 = (p1 - p0) ^ (p2 - p0);
        area = normal0.norm() * 0.5;

        Point3d axis0_1 = p1 - p0;
        axis0_1 /= axis0_1.norm();
        Point3d axis0_2 = normal0 ^ axis0_1;
        axis0_2 /= axis0_2.norm();

        Point2d p[3];
        p[0] = Point2d(p0 * (axis0_2), p0 * (axis0_1));
        p[1] = Point2d(p1 * (axis0_2), p1 * (axis0_1));
        p[2] = Point2d(p2 * (axis0_2), p2 * (axis0_1));

        Matrix2d P;
        P[0][0] = p[0][0] - p[2][0];
        P[1][0] = p[0][1] - p[2][1];
        P[0][1] = p[1][0] - p[2][0];
        P[1][1] = p[1][1] - p[2][1];

        const double Mdet = det(P);
        if (fabs(Mdet) > EPS)
        {
                invRestMat = inverse(P);
                return true;
        }
        return false;
}

// ----------------------------------------------------------------------------------------------
bool PBD::solve_FEMTriangleConstraint(
        const Point3d &p0, double invMass0,
        const Point3d &p1, double invMass1,
        const Point3d &p2, double invMass2,
        const double &area,
        const Matrix2d &invRestMat,
        const double youngsModulusX,
        const double youngsModulusY,
        const double youngsModulusShear,
        const double poissonRatioXY,
        const double poissonRatioYX,
        Point3d &corr0, Point3d &corr1, Point3d &corr2)
{
        // Orthotropic elasticity tensor
        Matrix3d C;
        C[0][0] = youngsModulusX / (1 - poissonRatioXY*poissonRatioYX);
        C[0][1] = youngsModulusX*poissonRatioYX / (static_cast<double>(1.0) - poissonRatioXY*poissonRatioYX);
        C[1][1] = youngsModulusY / (1 - poissonRatioXY*poissonRatioYX);
        C[1][0]= youngsModulusY*poissonRatioXY / (1 - poissonRatioXY*poissonRatioYX);
        C[2][2] = youngsModulusShear;

        // Determine \partial x/\partial m_i
        Matrix3x2d F;
        const Point3d p13 = p0 - p2;
        const Point3d p23 = p1 - p2;
        F[0][0] = p13[0] * invRestMat[0][0] + p23[0] * invRestMat[1][0];
        F[0][1] = p13[0] * invRestMat[0][1] + p23[0] * invRestMat[1][1];
        F[1][0] = p13[1] * invRestMat[0][0] + p23[1] * invRestMat[1][0];
        F[1][1] = p13[1] * invRestMat[0][1] + p23[1] * invRestMat[1][1];
        F[2][0] = p13[2] * invRestMat[0][0] + p23[2] * invRestMat[1][0];
        F[2][1] = p13[2] * invRestMat[0][1] + p23[2] * invRestMat[1][1];

        // epsilon = 0.5(F^T * F - I)
        Matrix2d epsilon;
        epsilon[0][0] = 0.5*(F[0][0] * F[0][0] + F[1][0] * F[1][0] + F[2][0] * F[2][0] - 1);		// xx
        epsilon[1][1] = 0.5*(F[0][1] * F[0][1] + F[1][1] * F[1][1] + F[2][1] * F[2][1] - 1);		// yy
        epsilon[0][1] = 0.5*(F[0][0] * F[0][1] + F[1][0] * F[1][1] + F[2][0] * F[2][1]);			// xy
        epsilon[1][0] = epsilon[0][1];

        // P(F) = det(F) * C*E * F^-T => E = green strain
        Matrix2d stress;
        stress[0][0] = C[0][0] * epsilon[0][0] + C[0][1] * epsilon[1][1] + C[0][2] * epsilon[0][1];
        stress[1][1] = C[1][0] * epsilon[0][0] + C[1][1] * epsilon[1][1] + C[1][2] * epsilon[0][1];
        stress[0][1] = C[2][0] * epsilon[0][0] + C[2][1] * epsilon[1][1] + C[2][2] * epsilon[0][1];
        stress[1][0] = stress[0][1];

        Matrix3x2d  piolaKirchhoffStres = F * stress;

        double psi = 0.0;
        for (unsigned char j = 0; j < 2; j++)
                for (unsigned char k = 0; k < 2; k++)
                        psi += epsilon(j,k) * stress(j,k);
        psi = 0.5*psi;
        double energy = area*psi;

        // compute gradient
        Matrix3x2d H = area * piolaKirchhoffStres * transpose(invRestMat);
        Point3d gradC[3];
        for (unsigned char j = 0; j < 3; ++j)
        {
                gradC[0][j] = H[j][0];
                gradC[1][j] = H[j][1];
        }
        gradC[2] = -gradC[0] - gradC[1];


        double sum_normGradC = invMass0 * gradC[0].norm() * gradC[0].norm() ;
        sum_normGradC += invMass1 *  gradC[1].norm() * gradC[1].norm();
        sum_normGradC += invMass2 *  gradC[2].norm() * gradC[2].norm();

        // exit early if required
        if (fabs(sum_normGradC) > EPS)
        {
                // compute scaling factor
                const double s = energy / sum_normGradC;

                // update positions
                corr0 = -(s*invMass0) * gradC[0];
                corr1 = -(s*invMass1) * gradC[1];
                corr2 = -(s*invMass2) * gradC[2];

                return true;
        }

        return false;
}


bool PBD::initialize(QWidget* parent) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("CellTissueProcess::initialize No current mesh"));

    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    CCStructure& cs =mesh->ccStructure("Tissue");

    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::FaceDataAttr &faceAttr =
      mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    Tissue::VertexDataAttr &vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");


    printStats = parm("Print Stats") == "True";
    velocityUpdate = parm("Velocity Update");
    maxVelocity = parm("Max Velocity").toDouble();
    staticDamping = parm("Velocity Static Damping").toDouble();
    viscousDamping = parm("Viscous Damping").toDouble();
    PBDiterations = parm("PBD iterations").toInt();
    PBDpressureStiffness = parm("PBD pressure stiffness").toDouble();
    PBDpressureAlpha =  parm("PBD pressure alpha").toDouble();
    PBDdistanceStiffness = parm("PBD distance stiffness").toDouble();
    PBDdistanceAlpha = parm("PBD distance alpha").toDouble();
    PBDbendingStiffness = parm("PBD bending stiffness").toDouble();
    PBDstrainStiffness = parm("PBD strain stiffness").toDouble();
    //PBDstrainAlpha = parm("PBD strain alpha").toDouble();
    PBDstrainShearStiffness = parm("PBD strain shear stiffness").toDouble();
    PBDstrainNormalize = parm("PBD strain normalize") == "True";
    PBDshapeStiffness = parm("PBD shape stiffness").toDouble();
    PBDshapeAllowStretching = parm("PBD shape allow stretching") == "True";

    if(printStats)
        mdxInfo << "PBD::initialize" << endl;


    // shape matching constrain initial
    int shape_init_n = 0;
    std::vector<int> labels;
    for(auto c : cellAttr) {
        Tissue::CellData& cD = cellAttr[c.first];
        if(!cD.shapeInit) {
            cD.restX0.clear();
            std::vector<double> invMasses;
            for(CCIndex v : cD.cellVertices) {
                cD.restX0.push_back(Point3d(indexAttr[v].pos));
                invMasses.push_back(vMAttr[v].invmass);
            }
            init_ShapeMatchingConstraint(cD.restX0, invMasses, cD.restX0.size(), cD.restCm, cD.invRestMat);
            for(CCIndex f : *cD.cellFaces) {
                Tissue::FaceData &fD = faceAttr[f];
                std::vector<CCIndex> vs = faceVertices(cs, f);
                fD.restPos[0] = indexAttr[vs[0]].pos;
                fD.restPos[1] = indexAttr[vs[1]].pos;
                fD.restPos[2] = indexAttr[vs[2]].pos;
            }
            shape_init_n++;
            labels.push_back(cD.label);
            cD.shapeInit = true;
        }
    }

    if(printStats) {
        mdxInfo << "PBD::initialize: Shape reinitialized for " << shape_init_n << " cells: ";
        for(int label : labels)
           mdxInfo << label << " " ;
        mdxInfo << endl;
    }


    return true;
}

// Rest pos attributes creates problem when starting from a saved mesh. Rewind the PBD process to fix it
bool PBD::rewind(QWidget* parent) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("PBD::initialize No current mesh"));

    CCStructure& cs = mesh->ccStructure("Tissue");
    CCIndexDataAttr& indexAttr = mesh->indexAttr();

    Tissue::FaceDataAttr &faceAttr =
      mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    for(CCIndex f : cs.faces()) {
        faceAttr[f].invRestMat = 0;
        for(int i = 0; i < 3; i++) {
            std::vector<CCIndex> vs = faceVertices(cs, f);
            faceAttr[f].restPos[i] = indexAttr[vs[i]].pos;
        }
    }
    return true;
}

// does not seems to work
void PBD::momentumPreservation(const CCStructure& cs, const CCIndexDataAttr &indexAttr, Tissue::VertexDataAttr &vMAttr, double K) {
  const CCIndexVec &vertices = cs.vertices();
  Point3d com, v_com, P_pbd, P_r, L_r, L_pbd, L;
  Matrix3d I;
  ////#pragma omp parallel for
  for (uint i = 0; i < vertices.size(); i++) {
    CCIndex v = vertices[i];
    com += indexAttr[v].pos;
    v_com += vMAttr[v].velocity;
  }
  com /= vertices.size();
  v_com /= vertices.size();
  ////#pragma omp parallel for
  for (uint i = 0; i < vertices.size(); i++) {
    CCIndex v = vertices[i];
    Point3d r = indexAttr[v].pos - com;
    Matrix3d r_star;
    r_star[0] = Point3d(0, -r[2], r[1]);
    r_star[1] = Point3d(r[2], 0, -r[0]);
    r_star[2] = Point3d(-r[1], r[0], 0);
    I += r_star * transpose(r_star);
    P_pbd += vMAttr[v].velocity;
    L_pbd += r ^ vMAttr[v].velocity;
    P_r += (vMAttr[v].velocity - vMAttr[v].prevVelocity);
    L_r += (indexAttr[v].pos - com) * vMAttr[v].velocity ;
  }
  Point3d omega_cor = inverse(I) * L_pbd;
  ////#pragma omp parallel for
  for (uint i = 0; i < vertices.size(); i++) {
    CCIndex v = vertices[i];
    Point3d delta_v =
        v_com + (omega_cor ^ (indexAttr[v].pos - com)) - vMAttr[v].velocity;
    vMAttr[v].velocity += delta_v * K;
  }
}

void PBD::dampVelocity(CCIndex v, const CCStructure& cs, const CCIndexDataAttr& indexAttr,
                                    Tissue::VertexDataAttr &vMAttr, Tissue::VertexData& vD) {
    Point3d damping = vD.velocity * staticDamping * Dt  ;
    Point3d vPos = indexAttr[v].pos;
    for(CCIndex vn : cs.neighbors(v)) {
        Tissue::VertexData& vDn = vMAttr[vn];
        Point3d nPos = indexAttr[vn].pos;
        damping += 0.5 * viscousDamping * Dt * ((vD.velocity - vDn.velocity) / norm(vPos - nPos));
    }
    vD.velocity -= damping;
    vD.dampedVelocity = damping;
    if(norm(vD.velocity) > maxVelocity)
        vD.velocity -= (norm(vD.velocity) - maxVelocity) * (vD.velocity / norm(vD.velocity));

}

// Calculate forces and velocities
void PBD::semiImplicitEuler(CCIndex v, const CCStructure& cs,
                                          CCIndexDataAttr& indexAttr, Tissue::VertexDataAttr &vMAttr, Tissue::VertexData& vD) {

    vD.prevVelocity = vD.velocity;
    vD.velocity += Dt * vD.force * vD.invmass;
    dampVelocity(v, cs, indexAttr, vMAttr, vD);
    vD.lastPos = vD.prevPos;
    vD.prevPos = indexAttr[v].pos;
    indexAttr[v].pos += vD.velocity * Dt;
    indexAttr[v].pos[2] = 0;
}




void PBD::solve() {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("CellTissueProcess::initialize No current mesh"));

    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    Tissue::CellDataAttr& cellAttr = mesh->attributes().attrMap<int, Tissue::CellData>("CellData");
    Tissue::FaceDataAttr &faceAttr =
      mesh->attributes().attrMap<CCIndex, Tissue::FaceData>("FaceData");
    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
    Tissue::EdgeDataAttr& edgeAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::EdgeData>("EdgeData");
    CCStructure& cs = mesh->ccStructure("Tissue");

    // Update proposal positions and initialize constrain corrections
    ////#pragma omp parallel for
    for(uint i = 0; i < cs.vertices().size(); i++) {
        CCIndex v = cs.vertices()[i];
        Tissue::VertexData& vD = vMAttr[v];
        semiImplicitEuler(v, cs, indexAttr, vMAttr, vD);
        vD.corrections.clear();
    }

    // Momentum preservation
    if(parm("Momentum Preservation K").toDouble() > 0)
        momentumPreservation(cs, indexAttr, vMAttr, parm("Momentum Preservation K").toDouble());

    // Apply PBD constrains
    double lambda_pressure[cellAttr.size()];
    for(uint i = 0; i < cellAttr.size(); i++)
        lambda_pressure[i] = 0;
    double lambda_distance[cs.edges().size()];
    for(uint i = 0; i < cs.edges().size(); i++)
        lambda_distance[i] = 0;
    for(int inter = 0; inter < PBDiterations; inter++) {
        // Pressure constraint
        if(PBDpressureStiffness > 0) {
            ////#pragma omp parallel for
            for(uint i = 0; i < cellAttr.size(); i++) {
                auto it = cellAttr.begin();
                advance(it, i);
                Tissue::CellData& cD = it->second;
                double compliance = PBDpressureAlpha / (Dt*Dt);
                std::map<CCIndex, Point3d> corr;
                if(!solve_PressureConstraint(cs, indexAttr, vMAttr, *cD.cellFaces, cD.restArea, cD.pressure, corr, lambda_pressure[i], compliance)) {
                    mdxInfo << "solvePBD: Failed solving pressure constraint for cell: " << cD.label << endl;
                    continue;
                }
                double stiffness = PBDpressureStiffness;
                if(stiffnessCorrection)
                    stiffness = 1 - pow(1 - stiffness, 1. / (inter + 1));
                for(CCIndex v : cD.cellVertices) {
                    Tissue::VertexData& vD = vMAttr[v];
                    Point3d corrv = /*1./vD.clusters **/ corr[v] * stiffness;
                    vD.corrections["pressure"] += corrv;
                    indexAttr[v].pos += corrv;
                }
            }
        }

        if(PBDdistanceStiffness > 0) {
            // Distance constraint
            ////#pragma omp parallel for shared(lockr, lockw)
            for(uint i = 0; i < cs.edges().size(); i++) {
                CCIndex e = cs.edges()[i];
                Tissue::EdgeData& eD = edgeAttr[e];
                auto eb = cs.edgeBounds(e);
                Point3d p0 = indexAttr[eb.first].pos;
                Point3d p1 = indexAttr[eb.second].pos;
                double m0 = vMAttr[eb.first].invmass;
                double m1 = vMAttr[eb.second].invmass;
                Point3d corr0, corr1;
                double stiffness = PBDdistanceStiffness;
                if(stiffnessCorrection)
                    stiffness = 1 - pow(1 - stiffness, 1. / (inter + 1));
                double compliance = PBDdistanceAlpha / (Dt*Dt);
                if(PBDdistanceAlpha == 0) {
                    if(!solve_DistanceConstraint(p0, m0, p1, m1, eD.restLength, eD.cStiffness, eD.eStiffness, corr0, corr1)) {
                        mdxInfo << "solvePBD: Failed solving distance constraint for edge: " << e << endl;
                        continue;
                    }
                } else {
                    if(!solve_DistanceConstraint_XPBD(p0, m0, p1, m1, eD.restLength, eD.cStiffness, eD.eStiffness, corr0, corr1, lambda_distance[i], compliance)){
                        mdxInfo << "solvePBD: Failed solving distance constraint for edge: " << e << endl;
                        continue;
                    }
                }
                corr0 *= stiffness;
                corr1 *= stiffness;
                vMAttr[eb.first].corrections["distance"] += corr0;
                vMAttr[eb.second].corrections["distance"] += corr1;
                indexAttr[eb.first].pos += corr0;
                indexAttr[eb.second].pos += corr1;
            }
        }

        if(PBDstrainStiffness > 0) {
            // Strains Constraint
            ////#pragma omp parallel for schedule(static)
            for(uint i = 0; i < cs.faces().size(); i++) {
                CCIndex f = cs.faces()[i];
                Tissue::FaceData& fD = faceAttr[f];
                std::vector<CCIndex> vs = faceVertices(cs, f);
                Point3d p0 = indexAttr[vs[0]].pos;
                Point3d p1 = indexAttr[vs[1]].pos;
                Point3d p2 = indexAttr[vs[2]].pos;
                double faceAngle = mdx::angle(Point3d(1,0,0), fD.a1/norm(fD.a1), Point3d(0,0,-1));
                Point3d x0 = Rotate(fD.restPos[0], faceAngle);
                Point3d x1 = Rotate(fD.restPos[1], faceAngle);
                Point3d x2 = Rotate(fD.restPos[2], faceAngle);
                if(!init_StrainTriangleConstraint(x0, x1, x2, fD.invRestMat)) {
                                    mdxInfo << "solvePBD: Failed initialization of strain constraint for face: " << f << endl;
                                    continue;
                }
                /*
                if(!init_FEMTriangleConstraint(x0, x1, x2, fD.restAreaFace, fD.invRestMat)) {
                    mdxInfo << "solvePBD: Failed initialization of strain constraint for face: " << f << endl;
                    continue;
                }
                */
                double m0 = vMAttr[vs[0]].invmass;
                double m1 = vMAttr[vs[1]].invmass;
                double m2 = vMAttr[vs[2]].invmass;
                //double compliance = PBDstrainAlpha / (Dt*Dt);
                Point3d corr0, corr1, corr2;
                double stiffnessXX = 0;
                double stiffnessYY = 0;
                double stiffnessXY = PBDstrainShearStiffness;
                stiffnessYY = norm(fD.a1) - norm(fD.a2);
                stiffnessXX = 0;//norm(fD.a2);
                if(!solve_StrainTriangleConstraint(p0, m0, p1, m1, p2, m2, fD.invRestMat, stiffnessXX, stiffnessYY, stiffnessXY,
                                                   PBDstrainNormalize, PBDstrainNormalize, corr0, corr1, corr2)) {
                    mdxInfo << "solvePBD: Failed solving strain constraint for face: " << f << endl;
                    continue;
                }
                /*
                if(!solve_FEMTriangleConstraint(p0, m0, p1, m1, p2, m2, fD.restAreaFace, fD.invRestMat, stiffnessXX, stiffnessYY, stiffnessXY,
0.01, 0.01,    corr0, corr1, corr2))   {
                    mdxInfo << "solvePBD: Failed solving strain constraint for face: " << f << endl;
                    continue;
                }
                */
                double stiffness = PBDstrainStiffness;
                if(stiffnessCorrection)
                    stiffness = 1 - pow(1 - stiffness, 1. / (inter + 1));
                corr0 *= stiffness;
                corr1 *= stiffness;
                corr2 *= stiffness;
                vMAttr[vs[0]].corrections["strain"] += corr0;
                vMAttr[vs[1]].corrections["strain"] += corr1;
                vMAttr[vs[2]].corrections["strain"] += corr2;
                indexAttr[vs[0]].pos += corr0;
                indexAttr[vs[1]].pos += corr1;
                indexAttr[vs[2]].pos += corr2;
            }
            //cout << avg_a1 / cs.faces().size() << " " << avg_a2 / cs.faces().size() << endl;

        }
        if(PBDshapeStiffness> 0) {
            // Shape matching constraint
            ////#pragma omp parallel for schedule(static)
            for(uint i = 0; i < cellAttr.size(); i++) {
                auto it = cellAttr.begin();
                advance(it, i);
                Tissue::CellData& cD = it->second;
                std::vector<Point3d> X;
                std::vector<double> invMasses;
                for(CCIndex v : cD.cellVertices) {
                    X.push_back(indexAttr[v].pos);
                    invMasses.push_back(vMAttr[v].invmass);
                }
                std::vector<Point3d> corr;
                Matrix3d rot;
                if(!solve_ShapeMatchingConstraint(cD.restX0,
                                              X,
                                              invMasses,
                                              cD.cellVertices.size(),
                                              cD.restCm,
                                              cD.invRestMat,
                                              PBDshapeAllowStretching,
                                              corr,
                                              &rot))
                {
                    mdxInfo << "solvePBD: Failed solving shape constraint for cell: " << cD.label << endl;
                    continue;
                }
                double stiffness = PBDshapeStiffness;
                if(stiffnessCorrection)
                    stiffness = 1 - pow(1 - stiffness, 1. / (inter + 1));
                int k = 0;
                for(CCIndex v : cD.cellVertices) {
                    Tissue::VertexData& vD = vMAttr[v];
                    Point3d corrv = /*1./vD.clusters * */ corr[k] * stiffness; // dividing by clusters does not preserve linear momentum
                    vD.corrections["shape"] += corrv;
                    indexAttr[v].pos += corrv;
                    k++;
                }
            }
        }
        if(PBDbendingStiffness > 0) {
            // Bending constraint
            ////#pragma omp parallel for
            for(uint i = 0; i < cellAttr.size(); i++) {
                auto it = cellAttr.begin();
                advance(it, i);
                Tissue::CellData& cD = it->second;
                int n = cD.perimeterVertices.size();
                for(int i = 0; i < n; i++) {
                    CCIndex v = cD.perimeterVertices[i];
                    CCIndex prev = (i>0) ? cD.perimeterVertices[i-1] : cD.perimeterVertices[n-1];
                    CCIndex next = (i<n-1) ? cD.perimeterVertices[i+1] : cD.perimeterVertices[0];
                    CCIndex e1, e2;
                    for(CCIndex e : cs.incidentCells(v, 1))
                        if(find(cD.perimeterEdges.begin(), cD.perimeterEdges.end(), e) != cD.perimeterEdges.end()){
                            if(e1.isPseudocell())
                                e1 = e;
                            else
                                e2 = e;
                        }
                    if(e1.isPseudocell() || e2.isPseudocell())
                        throw(QString("Somwthing wrong here"));
                    //double restAngle = fmod(cD.perimeterAngles[v]-M_PI, M_PI);
                    double restAngle = fmod(vMAttr[v].angle[make_pair(e1, e2)]-M_PI, M_PI);
                    if(restAngle < 0)
                        restAngle *= -1;
                    Point3d p0 = indexAttr[v].pos;
                    Point3d p1 = indexAttr[v].pos + Point3d(0,0,-1);
                    Point3d p2 = indexAttr[prev].pos;
                    Point3d p3 = indexAttr[next].pos;
                    double m0 = vMAttr[v].invmass;
                    double m1 = vMAttr[v].invmass;
                    double m2 = vMAttr[prev].invmass;
                    double m3 = vMAttr[next].invmass;
                    Point3d corr0, corr1, corr2, corr3;
                    double stiffness = PBDbendingStiffness;
                    if(stiffnessCorrection)
                                        stiffness = 1 - pow(1 - stiffness, 1. / (inter + 1));
                    solve_DihedralConstraint(p2, m2, p3, m3, p0, m0, p1, m1, restAngle, stiffness, corr2, corr3, corr0, corr1, indexAttr[v].selected);
                    vMAttr[v].corrections["bending"] += corr0;
                    vMAttr[prev].corrections["bending"] += corr2;
                    vMAttr[next].corrections["bending"] += corr3;
                    indexAttr[v].pos += corr0;
                    indexAttr[prev].pos += corr2;
                    indexAttr[next].pos += corr3;
                }
            }
        }

        //#pragma omp parallel for
        for(uint i = 0; i < cs.vertices().size(); i++) {
            CCIndex v = cs.vertices()[i];
            Tissue::VertexData& vD = vMAttr[v];
            indexAttr[v].pos[2] = 0;
            if(parm("Substrate Fixed") == "True" && (vD.substrate || vD.source)) {
                indexAttr[v].pos -= indexAttr[v].pos - vD.prevPos;
                vD.corrections["substrate/source"] -= indexAttr[v].pos - vD.prevPos;
            }
        }
    }

}



void PBD::update(double Dt) {
    Mesh* mesh = getMesh("Mesh 1");
    if(!mesh)
        throw(QString("CellTissueProcess::initialize No current mesh"));

    Tissue::VertexDataAttr& vMAttr =
        mesh->attributes().attrMap<CCIndex, Tissue::VertexData>("VertexData");
    CCIndexDataAttr& indexAttr = mesh->indexAttr();
    CCStructure& cs = mesh->ccStructure("Tissue");

    this->Dt = Dt;

    // call the PBD solver
    solve();

    mesh->updatePositions("Tissue");

    // update velocity
    //#pragma omp parallel for
    for(uint i = 0; i < cs.vertices().size(); i++) {
        CCIndex v = cs.vertices()[i];
        Tissue::VertexData& vD = vMAttr[v];
        vD.prevVelocity = vD.velocity;
        if(velocityUpdate == "First Order")
            vD.velocity = (1.0 / Dt) * (indexAttr[v].pos - vD.prevPos) ;
        else if (velocityUpdate == "Second Order")
            vD.velocity = (1.0 / Dt) * (1.5 * indexAttr[v].pos - 2.0 * vD.prevPos + 0.5 * vD.lastPos);
        else
            throw(QString("Wrong velocity order"));
    }

    // Get some debug info
    double avg_velocity = 0, avg_damping = 0;
    std::map<const char*, double> avg_corrections;
    Point3d linearMomentum, angularMomentum;
    for(uint i = 0; i < cs.vertices().size(); i++) {
        CCIndex v = cs.vertices()[i];
        Tissue::VertexData& vD = vMAttr[v];
        avg_velocity += norm(vD.velocity);
        avg_damping += norm(vD.dampedVelocity);
        for(auto deltap : vD.corrections) {
            const char* constraint = deltap.first;
            if(QString(constraint) == "substrate/source" || QString(constraint) == "substrate/source") continue;
            Point3d corr = deltap.second;
            avg_corrections[constraint] += norm(corr);
            linearMomentum += 1./vD.invmass * corr;
            angularMomentum += indexAttr[v].pos ^ (1./vD.invmass * corr);
        }
    }

    // Print some debug info
    if(printStats && debug_step > 10) {
        mdxInfo << "PDB: ";
        mdxInfo << " Avg velocity: " << avg_velocity / cs.vertices().size()
                << " Avg damping: " << avg_damping / cs.vertices().size();
        for(auto corrs : avg_corrections)
        mdxInfo << " Avg " << corrs.first << ": " << corrs.second / cs.vertices().size() << " ";
        mdxInfo << " Linear Momentum: " << norm(linearMomentum) << " ";
        mdxInfo << " Angular Momentum: " << norm(angularMomentum) << endl;
        debug_step = 0;
    }
    debug_step++;

}
