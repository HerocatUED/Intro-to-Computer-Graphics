#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <numeric>
#include <glm/gtx/vector_angle.hpp>
#include <spdlog/spdlog.h>

#include "Engine/Scene.h"
#include "Labs/FinalProject/Ray.h"

#define INF 1e+10f
#define PI 3.141592653589793f

namespace VCX::Labs::Rendering {

    constexpr float EPS1 = 1e-2f; // distance to prevent self-intersection
    constexpr float EPS2 = 1e-8f; // angle for parallel judgement
    constexpr float EPS3 = 1e-4f; // relative distance to enlarge kdtree
    constexpr float RussianRoulette = 0.8f;

    glm::vec4 GetTexture(Engine::Texture2D<Engine::Formats::RGBA8> const & texture, glm::vec2 const & uvCoord);

    glm::vec4 GetAlbedo(Engine::Material const & material, glm::vec2 const & uvCoord);

    glm::vec3 FR(float const & roughness, float const & shininess, glm::vec3 const & kd, glm::vec3 const & ks, glm::vec3 const & wi, glm::vec3 const & w0, glm::vec3 const & N);

    float smith(float const & roughness, glm::vec3 const & wi, glm::vec3 const & w0, glm::vec3 const & h);

    float smithG1(float const & roughness, glm::vec3 const & w, glm::vec3 const & h);

    float D(float const & roughness, glm::vec3 const & h, glm::vec3 const & N);

    float schlick(float const & shininess, glm::vec3 const & wi, glm::vec3 const & h);

    float fresnel(float const & shininess, glm::vec3 const & wi, glm::vec3 const & N);

    glm::vec3 SampleLight(float const & roughness, glm::vec3 const & ks, glm::vec3 const & wi, glm::vec3 const & N, float & pdf, bool const & enableImprotanceSampling);

    glm::vec3 SamplePoint(VCX::Engine::Light const & light);

    glm::vec3 toWorld(glm::vec3 const & a, glm::vec3 const & N);

    bool isLight(glm::vec3 const & pos, VCX::Engine::Light const & light);

    struct Intersection {
        float t, u, v; // ray parameter t, barycentric coordinates (u, v)
    };

    bool IntersectTriangle(Intersection & output, Ray const & ray, glm::vec3 const & p1, glm::vec3 const & p2, glm::vec3 const & p3);

    struct RayHit {
        bool              IntersectState;
        Engine::BlendMode IntersectMode;
        glm::vec3         IntersectPosition;
        glm::vec3         IntersectNormal;
        glm::vec4         IntersectAlbedo;   // [Albedo   (vec3), Alpha     (float)]
        glm::vec4         IntersectMetaSpec; // [Specular (vec3), Shininess (float)]
    };

    struct TrivialRayIntersector {
        Engine::Scene const * InternalScene = nullptr;

        TrivialRayIntersector() = default;

        void InitScene(Engine::Scene const * scene) {
            InternalScene = scene;
        }

        RayHit IntersectRay(Ray const & ray) const {
            RayHit result;
            if (! InternalScene) {
                spdlog::warn("VCX::Labs::Rendering::RayIntersector::IntersectRay(..): uninitialized intersector.");
                result.IntersectState = false;
                return result;
            }
            int          modelIdx, meshIdx;
            Intersection its;
            float        tmin     = 1e7, umin, vmin;
            int          maxmodel = InternalScene->Models.size();
            for (int i = 0; i < maxmodel; ++i) {
                auto const & model  = InternalScene->Models[i];
                int          maxidx = model.Mesh.Indices.size();
                for (int j = 0; j < maxidx; j += 3) {
                    std::uint32_t const * face = model.Mesh.Indices.data() + j;
                    glm::vec3 const &     p1   = model.Mesh.Positions[face[0]];
                    glm::vec3 const &     p2   = model.Mesh.Positions[face[1]];
                    glm::vec3 const &     p3   = model.Mesh.Positions[face[2]];
                    if (! IntersectTriangle(its, ray, p1, p2, p3)) continue;
                    if (its.t < EPS1 || its.t > tmin) continue;
                    tmin = its.t, umin = its.u, vmin = its.v, modelIdx = i, meshIdx = j;
                }
            }
            if (tmin == 1e7) {
                result.IntersectState = false;
                return result;
            }
            auto const &          model     = InternalScene->Models[modelIdx];
            auto const &          normals   = model.Mesh.IsNormalAvailable() ? model.Mesh.Normals : model.Mesh.ComputeNormals();
            auto const &          texcoords = model.Mesh.IsTexCoordAvailable() ? model.Mesh.TexCoords : model.Mesh.GetEmptyTexCoords();
            std::uint32_t const * face      = model.Mesh.Indices.data() + meshIdx;
            glm::vec3 const &     p1        = model.Mesh.Positions[face[0]];
            glm::vec3 const &     p2        = model.Mesh.Positions[face[1]];
            glm::vec3 const &     p3        = model.Mesh.Positions[face[2]];
            glm::vec3 const &     n1        = normals[face[0]];
            glm::vec3 const &     n2        = normals[face[1]];
            glm::vec3 const &     n3        = normals[face[2]];
            glm::vec2 const &     uv1       = texcoords[face[0]];
            glm::vec2 const &     uv2       = texcoords[face[1]];
            glm::vec2 const &     uv3       = texcoords[face[2]];
            result.IntersectState           = true;
            auto const & material           = InternalScene->Materials[model.MaterialIndex];
            result.IntersectMode            = material.Blend;
            result.IntersectPosition        = (1.0f - umin - vmin) * p1 + umin * p2 + vmin * p3;
            result.IntersectNormal          = (1.0f - umin - vmin) * n1 + umin * n2 + vmin * n3;
            glm::vec2 uvCoord               = (1.0f - umin - vmin) * uv1 + umin * uv2 + vmin * uv3;
            result.IntersectAlbedo          = GetAlbedo(material, uvCoord);
            result.IntersectMetaSpec        = GetTexture(material.MetaSpec, uvCoord);

            return result;
        }
    };

    /* Optional: write your own accelerated intersector here */

    typedef struct Triangle {
        glm::vec3 center;
        int       modelIdx, meshIdx;
        Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c, int x, int y) {
            modelIdx = x, meshIdx = y;
            center = (a + b + c) / glm::vec3(3, 3, 3);
        }
    } Triangle;

    bool cmpx(const Triangle & t1, const Triangle & t2);
    bool cmpy(const Triangle & t1, const Triangle & t2);
    bool cmpz(const Triangle & t1, const Triangle & t2);

    struct BVHNode {
        BVHNode * left  = NULL;
        BVHNode * right = NULL;
        int       n = 0, index = -1; // 叶子节点信息
        glm::vec3 AA = glm::vec3(0.0f), BB = glm::vec3(0.0f); // 碰撞盒
    };

    struct HitResult {
        int          idx = -1;
        Intersection its;
    };

    struct BVHRayIntersector {
        Engine::Scene const * InternalScene = nullptr;
        std::vector<Triangle> triangles;
        BVHNode *             root;

        BVHRayIntersector() = default;

        void my_destroy() {
            destroy(root);
            triangles.clear();
        }

        void destroy(BVHNode * n) {
            if (n->left)
                destroy(n->left);
            if (n->right)
                destroy(n->right);
            delete n;
            return;
        }

        BVHNode * buildBVHwithSAH(std::vector<Triangle> & triangles, int l, int r, int n) {
            BVHNode * node = new BVHNode();
            node->AA       = glm::vec3(INF, INF, INF);
            node->BB       = glm::vec3(-INF, -INF, -INF);
            for (int i = l; i <= r; i++) {
                // 最小点 AA
                auto const &          model = InternalScene->Models[triangles[i].modelIdx];
                std::uint32_t const * face  = model.Mesh.Indices.data() + triangles[i].meshIdx;
                glm::vec3 const &     p1    = model.Mesh.Positions[face[0]];
                glm::vec3 const &     p2    = model.Mesh.Positions[face[1]];
                glm::vec3 const &     p3    = model.Mesh.Positions[face[2]];
                float                 minx  = std::min(p1.x, std::min(p2.x, p3.x));
                float                 miny  = std::min(p1.y, std::min(p2.y, p3.y));
                float                 minz  = std::min(p1.z, std::min(p2.z, p3.z));
                node->AA.x                  = std::min(node->AA.x, minx);
                node->AA.y                  = std::min(node->AA.y, miny);
                node->AA.z                  = std::min(node->AA.z, minz);
                // 最大点 BB
                float maxx = std::max(p1.x, std::max(p2.x, p3.x));
                float maxy = std::max(p1.y, std::max(p2.y, p3.y));
                float maxz = std::max(p1.z, std::max(p2.z, p3.z));
                node->BB.x = std::max(node->BB.x, maxx);
                node->BB.y = std::max(node->BB.y, maxy);
                node->BB.z = std::max(node->BB.z, maxz);
            }
            // 不多于 n 个三角形 返回叶子节点
            if ((r - l + 1) <= n) {
                node->n     = r - l + 1;
                node->index = l;
                return node;
            }
            // 否则递归建树
            float Cost  = INF;
            int   Axis  = 0;
            int   Split = (l + r) / 2;
            for (int axis = 0; axis < 3; axis++) {
                if (axis == 0) std::sort(&triangles[l], &triangles[r] + 1, cmpx);
                if (axis == 1) std::sort(&triangles[l], &triangles[r] + 1, cmpy);
                if (axis == 2) std::sort(&triangles[l], &triangles[r] + 1, cmpz);
                std::vector<glm::vec3> leftMax(r - l + 1, glm::vec3(-INF, -INF, -INF));
                std::vector<glm::vec3> leftMin(r - l + 1, glm::vec3(INF, INF, INF));
                for (int i = l; i <= r; i++) {
                    auto const &          model = InternalScene->Models[triangles[i].modelIdx];
                    std::uint32_t const * face  = model.Mesh.Indices.data() + triangles[i].meshIdx;
                    glm::vec3 const &     p1    = model.Mesh.Positions[face[0]];
                    glm::vec3 const &     p2    = model.Mesh.Positions[face[1]];
                    glm::vec3 const &     p3    = model.Mesh.Positions[face[2]];
                    int                   bias  = (i == l) ? 0 : 1; // 第一个元素特殊处理
                    leftMax[i - l].x            = std::max(leftMax[i - l - bias].x, std::max(p1.x, std::max(p2.x, p3.x)));
                    leftMax[i - l].y            = std::max(leftMax[i - l - bias].y, std::max(p1.y, std::max(p2.y, p3.y)));
                    leftMax[i - l].z            = std::max(leftMax[i - l - bias].z, std::max(p1.z, std::max(p2.z, p3.z)));

                    leftMin[i - l].x = std::min(leftMin[i - l - bias].x, std::min(p1.x, std::min(p2.x, p3.x)));
                    leftMin[i - l].y = std::min(leftMin[i - l - bias].y, std::min(p1.y, std::min(p2.y, p3.y)));
                    leftMin[i - l].z = std::min(leftMin[i - l - bias].z, std::min(p1.z, std::min(p2.z, p3.z)));
                }
                std::vector<glm::vec3> rightMax(r - l + 1, glm::vec3(-INF, -INF, -INF));
                std::vector<glm::vec3> rightMin(r - l + 1, glm::vec3(INF, INF, INF));
                for (int i = r; i >= l; i--) {
                    auto const &          model = InternalScene->Models[triangles[i].modelIdx];
                    std::uint32_t const * face  = model.Mesh.Indices.data() + triangles[i].meshIdx;
                    glm::vec3 const &     p1    = model.Mesh.Positions[face[0]];
                    glm::vec3 const &     p2    = model.Mesh.Positions[face[1]];
                    glm::vec3 const &     p3    = model.Mesh.Positions[face[2]];
                    int                   bias  = (i == r) ? 0 : 1; // 第一个元素特殊处理
                    rightMax[i - l].x           = std::max(rightMax[i - l + bias].x, std::max(p1.x, std::max(p2.x, p3.x)));
                    rightMax[i - l].y           = std::max(rightMax[i - l + bias].y, std::max(p1.y, std::max(p2.y, p3.y)));
                    rightMax[i - l].z           = std::max(rightMax[i - l + bias].z, std::max(p1.z, std::max(p2.z, p3.z)));

                    rightMin[i - l].x = std::min(rightMin[i - l + bias].x, std::min(p1.x, std::min(p2.x, p3.x)));
                    rightMin[i - l].y = std::min(rightMin[i - l + bias].y, std::min(p1.y, std::min(p2.y, p3.y)));
                    rightMin[i - l].z = std::min(rightMin[i - l + bias].z, std::min(p1.z, std::min(p2.z, p3.z)));
                }
                // 遍历寻找分割
                float cost  = INF;
                int   split = l;
                for (int i = l ; i <= r - 1; i++) {
                    float lenx, leny, lenz;
                    // 左侧 [l, i]
                    glm::vec3 leftAA = leftMin[i - l];
                    glm::vec3 leftBB = leftMax[i - l];
                    lenx             = leftBB.x - leftAA.x;
                    leny             = leftBB.y - leftAA.y;
                    lenz             = leftBB.z - leftAA.z;
                    float leftS      = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
                    float leftCost   = leftS * (i - l + 1);
                    // 右侧 [i+1, r]
                    glm::vec3 rightAA = rightMin[i + 1 - l];
                    glm::vec3 rightBB = rightMax[i + 1 - l];
                    lenx              = rightBB.x - rightAA.x;
                    leny              = rightBB.y - rightAA.y;
                    lenz              = rightBB.z - rightAA.z;
                    float rightS      = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
                    float rightCost   = rightS * (r - i);
                    // 记录每个分割的最小答案
                    float totalCost = leftCost + rightCost;
                    if (totalCost < cost) {
                        cost  = totalCost;
                        split = i;
                    }
                }
                // 记录每个轴的最佳答案
                if (cost < Cost) {
                    Cost  = cost;
                    Axis  = axis;
                    Split = split;
                }
            }
            // 按最佳轴分割
            if (Axis == 0) std::sort(&triangles[l], &triangles[r] + 1, cmpx);
            if (Axis == 1) std::sort(&triangles[l], &triangles[r] + 1, cmpy);
            if (Axis == 2) std::sort(&triangles[l], &triangles[r] + 1, cmpz);
            node->left  = buildBVHwithSAH(triangles, l, Split, n);
            node->right = buildBVHwithSAH(triangles, Split + 1, r, n);
            return node;
        }

        void InitScene(Engine::Scene const * scene) {
            if (InternalScene) my_destroy();
            InternalScene = scene;
            int maxmodel  = InternalScene->Models.size();
            for (int i = 0; i < maxmodel; ++i) {
                auto const & model  = InternalScene->Models[i];
                int          maxidx = model.Mesh.Indices.size();
                for (int j = 0; j < maxidx; j += 3) {
                    std::uint32_t const * face = model.Mesh.Indices.data() + j;
                    glm::vec3 const &     p1   = model.Mesh.Positions[face[0]];
                    glm::vec3 const &     p2   = model.Mesh.Positions[face[1]];
                    glm::vec3 const &     p3   = model.Mesh.Positions[face[2]];
                    triangles.push_back(Triangle(p1, p2, p3, i, j));
                }
            }
            root = buildBVHwithSAH(triangles, 0, triangles.size() - 1, 8);
        }

        HitResult hitTriangleArray(Ray const & ray, std::vector<Triangle> const & triangles, int l, int r) const {
            HitResult    ans;
            Intersection tmp;
            bool         flag = false;
            for (int i = l; i <= r; i++) {
                auto const &          model = InternalScene->Models[triangles[i].modelIdx];
                std::uint32_t const * face  = model.Mesh.Indices.data() + triangles[i].meshIdx;
                glm::vec3 const &     p1    = model.Mesh.Positions[face[0]];
                glm::vec3 const &     p2    = model.Mesh.Positions[face[1]];
                glm::vec3 const &     p3    = model.Mesh.Positions[face[2]];
                if (IntersectTriangle(tmp, ray, p1, p2, p3) && tmp.t > EPS1) {
                    if (! flag) {
                        flag    = true;
                        ans.idx = i;
                        ans.its = tmp;
                    } else if (tmp.t < ans.its.t) {
                        ans.idx = i;
                        ans.its = tmp;
                    }
                }
            }
            return ans;
        }

        float hitAABB(Ray const & ray, glm::vec3 AA, glm::vec3 BB) const {
            glm::vec3 invdir = 1.0f / ray.Direction;
            glm::vec3 in     = (BB - ray.Origin) * invdir;
            glm::vec3 out    = (AA - ray.Origin) * invdir;
            glm::vec3 tmax   = max(in, out);
            glm::vec3 tmin   = min(in, out);
            float     t1     = std::min(tmax.x, std::min(tmax.y, tmax.z));
            float     t0     = std::max(tmin.x, std::max(tmin.y, tmin.z));
            return (t1 >= t0) ? ((t0 > 0.0) ? t0 : t1) : -1;
        }

        HitResult hitBVH(Ray const & ray, std::vector<Triangle> const & triangles, BVHNode const * root) const {
            if (root == NULL) return HitResult();
            if (root->n > 0) return hitTriangleArray(ray, triangles, root->index, root->n + root->index - 1);
            float d1 = -1, d2 = -1;
            if (root->left) d1 = hitAABB(ray, root->left->AA, root->left->BB);
            if (root->right) d2 = hitAABB(ray, root->right->AA, root->right->BB);
            HitResult r1, r2;
            if (d1 > 0) r1 = hitBVH(ray, triangles, root->left);
            if (d2 > 0) r2 = hitBVH(ray, triangles, root->right);
            if (r1.idx != -1 && r2.idx != -1) return r1.its.t < r2.its.t ? r1 : r2;
            else if (r1.idx != -1) return r1;
            else return r2;
        }

        RayHit IntersectRay(Ray const & ray) const {
            RayHit result;
            if (! InternalScene) {
                spdlog::warn("VCX::Labs::Rendering::RayIntersector::IntersectRay(..): uninitialized intersector.");
                result.IntersectState = false;
                return result;
            }
            HitResult res = hitBVH(ray, triangles, root);
            if (res.idx == -1) {
                result.IntersectState = false;
                return result;
            }
            int                   modelIdx = triangles[res.idx].modelIdx, meshIdx = triangles[res.idx].meshIdx;
            float                 umin = res.its.u, vmin = res.its.v;
            auto const &          model     = InternalScene->Models[modelIdx];
            auto const &          normals   = model.Mesh.IsNormalAvailable() ? model.Mesh.Normals : model.Mesh.ComputeNormals();
            auto const &          texcoords = model.Mesh.IsTexCoordAvailable() ? model.Mesh.TexCoords : model.Mesh.GetEmptyTexCoords();
            std::uint32_t const * face      = model.Mesh.Indices.data() + meshIdx;
            glm::vec3 const &     p1        = model.Mesh.Positions[face[0]];
            glm::vec3 const &     p2        = model.Mesh.Positions[face[1]];
            glm::vec3 const &     p3        = model.Mesh.Positions[face[2]];
            glm::vec3 const &     n1        = normals[face[0]];
            glm::vec3 const &     n2        = normals[face[1]];
            glm::vec3 const &     n3        = normals[face[2]];
            glm::vec2 const &     uv1       = texcoords[face[0]];
            glm::vec2 const &     uv2       = texcoords[face[1]];
            glm::vec2 const &     uv3       = texcoords[face[2]];
            result.IntersectState           = true;
            auto const & material           = InternalScene->Materials[model.MaterialIndex];
            result.IntersectMode            = material.Blend;
            result.IntersectPosition        = (1.0f - umin - vmin) * p1 + umin * p2 + vmin * p3;
            result.IntersectNormal          = (1.0f - umin - vmin) * n1 + umin * n2 + vmin * n3;
            glm::vec2 uvCoord               = (1.0f - umin - vmin) * uv1 + umin * uv2 + vmin * uv3;
            result.IntersectAlbedo          = GetAlbedo(material, uvCoord);
            result.IntersectMetaSpec        = GetTexture(material.MetaSpec, uvCoord);

            return result;
        }
    };

   //using RayIntersector = TrivialRayIntersector;
   using RayIntersector = BVHRayIntersector;

    glm::vec3 PathTrace(const RayIntersector & intersector, Ray ray, bool const & enableImprotanceSampling);

} // namespace VCX::Labs::Rendering
