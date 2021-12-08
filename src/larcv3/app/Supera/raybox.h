// This code is **based on** from https://www.scratchapixel.com/code.php?id=10&origin=/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes&src=1
// ... which is made available under GNU license

#ifndef LARCV_AABOX_H
#define LARCV_AABOX_H

#include "geometry.h" 
#include <cstdlib> 
#include "BBox.h"
//#include <random>

#if __has_include("larcv3/core/dataformat/Particle.h")
#define larcv larcv3
#include "larcv3/core/dataformat/ImageMeta.h"
#endif

    namespace larcv {

    template <typename T>
    class Ray 
    { 
    public: 
        Ray(const Vec3<T> &orig, const Vec3<T> &dir) : orig(orig), dir(dir)
        { 
            invdir = 1 / dir; 
            sign[0] = (invdir.x < 0); 
            sign[1] = (invdir.y < 0); 
            sign[2] = (invdir.z < 0); 
        } 
        Vec3<T> orig, dir; // ray orig and dir
        Vec3<T> invdir;
        int sign[3]; 
    }; 

    /// Axis-aligned bounding box.
    /// See https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
    template <typename T>
    class AABBox 
    { 
    public: 
        AABBox(const Vec3<T> &b0, const Vec3<T> &b1) { bounds[0] = b0, bounds[1] = b1; }
        AABBox(const T x0, const T y0, const T z0, const T x1, const T y1, const T z1)
        {
            bounds[0].x=x0; bounds[0].y=y0; bounds[0].z=z0;
            bounds[1].x=x1; bounds[1].y=y1; bounds[1].z=z1;
        }
        AABBox(const supera::BBox3D& box)
        {
            bounds[0].x=box.min_x(); bounds[0].y=box.min_y(); bounds[0].z=box.min_z(); 
            bounds[1].x=box.max_x(); bounds[1].y=box.max_y(); bounds[1].z=box.max_z();
        }
        #if __has_include("larcv3/core/dataformat/Particle.h")
        AABBox(const larcv3::ImageMeta<3>& box)
        {
            bounds[0].x=box.min(0); bounds[0].y=box.min(1); bounds[0].z=box.min(2); 
            bounds[1].x=box.max(0); bounds[1].y=box.max(1); bounds[1].z=box.max(2);
        }
        #endif

        /// Find locations of intersections of a given ray with this box.
        ///
        /// \param r   The ray
        /// \param t0  The intersection closer to the ray's origin
        /// \param t1  The intersection further from the ray's origin
        /// \return    The number of intersections:
        ///              If 0, t0 and t1 are meaningless
        ///              If 1, t0 is meaningless; t1 indicates the intersection point
        ///              If 2, both t0 and t1 are meaningful.
        int intersect(const Ray<T> &r, T& t0, T& t1) const
        { 
            T tmin, tmax, tymin, tymax, tzmin, tzmax;

            tmin = (bounds[r.sign[0]].x - r.orig.x) * r.invdir.x; 
            tmax = (bounds[1-r.sign[0]].x - r.orig.x) * r.invdir.x; 
            tymin = (bounds[r.sign[1]].y - r.orig.y) * r.invdir.y; 
            tymax = (bounds[1-r.sign[1]].y - r.orig.y) * r.invdir.y; 

            if ((tmin > tymax) || (tymin > tmax)) 
                return 0; 

            if (tymin > tmin) 
                tmin = tymin; 
            if (tymax < tmax) 
                tmax = tymax; 

            tzmin = (bounds[r.sign[2]].z - r.orig.z) * r.invdir.z; 
            tzmax = (bounds[1-r.sign[2]].z - r.orig.z) * r.invdir.z; 

            if ((tmin > tzmax) || (tzmin > tmax)) 
                return 0; 

            if (tzmin > tmin) 
                tmin = tzmin; 
            if (tzmax < tmax) 
                tmax = tzmax; 

            t0 = tmin; 
            t1 = tmax;

            if(t1 < 0) return 0;
            if(t0 < 0) return 1;
            return 2;
        } 

        bool contain(const Vec3<T>& pt) const {
            return (bounds[0].x <= pt.x && pt.x < bounds[1].x &&
                bounds[0].y <= pt.y && pt.y < bounds[1].y &&
                bounds[0].z <= pt.z && pt.z < bounds[1].z);
        }

        Vec3<T> bounds[2];
    }; 

}

#endif