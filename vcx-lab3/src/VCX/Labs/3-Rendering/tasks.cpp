#include "Labs/3-Rendering/tasks.h"

namespace VCX::Labs::Rendering {

    glm::vec4 GetTexture(Engine::Texture2D<Engine::Formats::RGBA8> const & texture, glm::vec2 const & uvCoord) {
        if (texture.GetSizeX() == 1 || texture.GetSizeY() == 1) return texture.At(0, 0);
        glm::vec2 uv      = glm::fract(uvCoord);
        uv.x              = uv.x * texture.GetSizeX() - .5f;
        uv.y              = uv.y * texture.GetSizeY() - .5f;
        std::size_t xmin  = std::size_t(glm::floor(uv.x) + texture.GetSizeX()) % texture.GetSizeX();
        std::size_t ymin  = std::size_t(glm::floor(uv.y) + texture.GetSizeY()) % texture.GetSizeY();
        std::size_t xmax  = (xmin + 1) % texture.GetSizeX();
        std::size_t ymax  = (ymin + 1) % texture.GetSizeY();
        float       xfrac = glm::fract(uv.x), yfrac = glm::fract(uv.y);
        return glm::mix(glm::mix(texture.At(xmin, ymin), texture.At(xmin, ymax), yfrac), glm::mix(texture.At(xmax, ymin), texture.At(xmax, ymax), yfrac), xfrac);
    }

    glm::vec4 GetAlbedo(Engine::Material const & material, glm::vec2 const & uvCoord) {
        glm::vec4 albedo       = GetTexture(material.Albedo, uvCoord);
        glm::vec3 diffuseColor = albedo;
        return glm::vec4(glm::pow(diffuseColor, glm::vec3(2.2)), albedo.w);
    }

    bool cmpx(const Triangle & t1, const Triangle & t2) {
        return t1.center.x < t2.center.x;
    }
    bool cmpy(const Triangle & t1, const Triangle & t2) {
        return t1.center.y < t2.center.y;
    }
    bool cmpz(const Triangle & t1, const Triangle & t2) {
        return t1.center.z < t2.center.z;
    }

    /******************* 1. Ray-triangle intersection *****************/
    bool IntersectTriangle(Intersection & output, Ray const & ray, glm::vec3 const & p1, glm::vec3 const & p2, glm::vec3 const & p3) {
        // your code here
        glm::vec3 edge1 = p2 - p1, edge2 = p3 - p1;
        glm::vec3 d    = glm::normalize(ray.Direction);
        glm::vec3 pvec = glm::cross(d, edge2);
        float     det  = glm::dot(edge1, pvec);
        if (det < EPS2) return false;
        glm::vec3 tvec = ray.Origin - p1;
        output.u       = glm::dot(tvec, pvec);
        if (output.u < 0.0f || output.u > det) return false;
        glm::vec3 qvec = glm::cross(tvec, edge1);
        output.v       = glm::dot(d, qvec);
        if (output.v < 0.0f || output.u + output.v > det) return false;
        output.t      = glm::dot(edge2, qvec);
        float inv_det = 1.0f / det;
        output.t *= inv_det;
        output.u *= inv_det;
        output.v *= inv_det;
    }

    glm::vec3 RayTrace(const RayIntersector & intersector, Ray ray, int maxDepth, bool enableShadow) {
        glm::vec3 color(0.0f);
        glm::vec3 weight(1.0f);

        for (int depth = 0; depth < maxDepth; depth++) {
            auto rayHit = intersector.IntersectRay(ray);
            if (! rayHit.IntersectState) return color;
            const glm::vec3 pos       = rayHit.IntersectPosition;
            const glm::vec3 n         = rayHit.IntersectNormal;
            const glm::vec3 kd        = rayHit.IntersectAlbedo;
            const glm::vec3 ks        = rayHit.IntersectMetaSpec;
            const float     alpha     = rayHit.IntersectAlbedo.w;
            const float     shininess = rayHit.IntersectMetaSpec.w * 256;

            glm::vec3 result(0.0f);
            /******************* 2. Whitted-style ray tracing *****************/
            // your code here
            result += 0.05f * kd;
            for (const Engine::Light & light : intersector.InternalScene->Lights) {
                glm::vec3 l;
                float     attenuation;
                /******************* 3. Shadow ray *****************/
                if (light.Type == Engine::LightType::Point) {
                    l           = light.Position - pos;
                    attenuation = 1.0f / glm::dot(l, l);
                    if (enableShadow) {
                        // your code here
                        auto  hit = intersector.IntersectRay(Ray(light.Position, -l));
                        float dis = glm::distance(pos, hit.IntersectPosition);
                        while (hit.IntersectState && dis > EPS1) {
                            if (hit.IntersectAlbedo.w >= 0.2f) break;
                            hit = intersector.IntersectRay(Ray(hit.IntersectPosition, -l));
                            dis = glm::distance(pos, hit.IntersectPosition);
                        }
                        if (hit.IntersectState && dis >= EPS1) continue;
                    }
                } else if (light.Type == Engine::LightType::Directional) {
                    l           = light.Direction;
                    attenuation = 1.0f;
                    if (enableShadow) {
                        // your code here
                        auto  hit = intersector.IntersectRay(Ray(pos + glm::normalize(l) * (1e+3f), -l));
                        float dis = glm::distance(pos, hit.IntersectPosition);
                        while (hit.IntersectState && dis > EPS1) {
                            if (hit.IntersectAlbedo.w >= 0.2f) break;
                            hit = intersector.IntersectRay(Ray(hit.IntersectPosition, -l));
                            dis = glm::distance(pos, hit.IntersectPosition);
                        }
                        if (hit.IntersectState && dis >= EPS1) continue;
                    }
                }

                /******************* 2. Whitted-style ray tracing *****************/
                // your code here
                l                  = glm::normalize(l);
                glm::vec3 d        = glm::normalize(-ray.Direction);
                glm::vec3 diffuse  = light.Intensity * kd * std::max(glm::dot(n, l), 0.0f);
                glm::vec3 specular = light.Intensity * ks * std::pow(std::max(glm::dot(n, glm::normalize(l + d)), 0.0f), shininess);
                result += attenuation * (diffuse + specular);
            }

            if (alpha < 0.9) {
                // refraction
                // accumulate color
                glm::vec3 R = alpha * glm::vec3(1.0f);
                color += weight * R * result;
                weight *= glm::vec3(1.0f) - R;

                // generate new ray
                ray = Ray(pos, ray.Direction);
            } else {
                // reflection
                // accumulate color
                glm::vec3 R = ks * glm::vec3(0.5f);
                color += weight * (glm::vec3(1.0f) - R) * result;
                weight *= R;

                // generate new ray
                glm::vec3 out_dir = ray.Direction - glm::vec3(2.0f) * n * glm::dot(n, ray.Direction);
                ray               = Ray(pos, out_dir);
            }
        }

        return color;
    }
} // namespace VCX::Labs::Rendering