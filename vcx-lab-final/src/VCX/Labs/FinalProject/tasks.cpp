#include "Labs/FinalProject/tasks.h"

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

    glm::vec3 SampleLight(float const & roughness, glm::vec3 const & ks, glm::vec3 const & wi, glm::vec3 const & N, float & pdf, bool const & enableImprotanceSampling) {
        static std::random_device                    dev;
        static std::mt19937                          rng(dev());
        static std::uniform_real_distribution<float> rand1(0.0f, 1.0f), rand2(0.0f, 1.0f);
        if (enableImprotanceSampling) {
            if (glm::dot(ks, ks) < EPS1) {
                // cosine-weighted for diffuse
                float     x_1 = rand1(rng), x_2 = rand2(rng);
                float     x = std::cos(2 * PI * x_2) * std::sqrt(x_1);
                float     y = std::sin(2 * PI * x_2) * std::sqrt(x_1);
                float     z = std::sqrt(1 - x_1);
                glm::vec3 localRay(x, y, z);
                pdf = z * std::sqrt(1 - z * z) / PI;
                return toWorld(localRay, N);
            } else {
                // ggx-BRDF for specular
                float     x_1 = rand1(rng), x_2 = rand2(rng);
                float     roughness2 = roughness * roughness;
                float     cos_theta2 = (1 - x_1) / (1 + x_1 * (roughness2 - 1));
                float     sin_theta  = std::sqrt(1 - cos_theta2);
                float     phi        = 2 * PI * x_2;
                float     x          = sin_theta * std::cos(phi);
                float     y          = sin_theta * std::sin(phi);
                float     z          = std::sqrt(cos_theta2);
                float     t          = (roughness2 - 1) * cos_theta2 + 1;
                float     pdf_h      = roughness2 * z * sin_theta / (PI * t * t);
                glm::vec3 localRay(x, y, z);
                glm::vec3 w0 = toWorld(localRay, N);
                glm::vec3 h  = glm::normalize(w0 + wi);
                pdf          = pdf_h / (4 * glm::dot(wi, h));
                return w0;
            }
        } else {
            // uniform sampling
            float     x_1 = rand1(rng), x_2 = rand2(rng);
            float     z = std::fabs(1.0f - 2.0f * x_1);
            float     r = std::sqrt(1.0f - z * z), phi = 2 * PI * x_2;
            glm::vec3 localRay(r * std::cos(phi), r * std::sin(phi), z);
            glm::vec3 w0 = toWorld(localRay, N);
            pdf          = glm::dot(w0, N) > 0 ? 0.5f / PI : 0.0f;
            return w0;
        }
    }

    glm::vec3 toWorld(glm::vec3 const & a, glm::vec3 const & N) {
        glm::vec3 B, C;
        if (std::fabs(N.x) > std::fabs(N.y)) {
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C            = glm::vec3(N.z * invLen, 0.0f, -N.x * invLen);
        } else {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C            = glm::vec3(0.0f, N.z * invLen, -N.y * invLen);
        }
        B = glm::cross(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

    bool isLight(glm::vec3 const & pos, VCX::Engine::Light const & light) {
        if (std::abs(pos.y - light.Position.y) < 1 && std::abs(pos.x - light.Position.x) < 65 && std::abs(pos.z - light.Position.z) < 57.5)
            return true;
        return false;
    }

    float fresnel(float const & shininess, glm::vec3 const & wi, glm::vec3 const & N) {
        float cosi = glm::dot(wi, N);
        float etai = 1, etat = shininess;
        if (cosi > 0) { std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(1 - cosi * cosi);
        // Total internal reflection
        if (sint >= 1) {
            return 1;
        } else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi       = fabsf(cosi);
            float Rs   = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp   = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            return (Rs * Rs + Rp * Rp) / 2;
        }
    }

    float schlick(float const & shininess, glm::vec3 const & wi, glm::vec3 const & N) {
        float cos_theta = glm::dot(wi, N);
        float etai = 1, etat = shininess;
        if (cos_theta > 0) { std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sin_theta = etai / etat * sqrtf(1 - cos_theta * cos_theta);
        if (sin_theta >= 1)
            return 1;
        float R0 = (1 - shininess) / (1 + shininess);
        R0       = R0 * R0;
        return R0 + (1 - R0) * std::pow((1 - cos_theta), 5);
    }

    float D(float const & roughness, glm::vec3 const & h, glm::vec3 const & N) {
        float cos_theta   = std::max(glm::dot(N, h), 0.0f);
        float cos_theta_2 = cos_theta * cos_theta;
        float roughness2  = roughness * roughness;
        // beckmann
        // float     exponential = (cos_theta_2 - 1) / (roughness2 * cos_theta_2);
        // float     d           = std::exp(exponential) / (PI * roughness2 * cos_theta_2 * cos_theta_2);
        // ggx
        float t = cos_theta_2 * (roughness2 - 1) + 1;
        float d = roughness2 / (PI * t * t);
        return d;
    }

    float smithG1(float const & roughness, glm::vec3 const & w, glm::vec3 const & h) {
        // Schlick-GGX
        float alpha     = roughness + 1;
        float k         = alpha * alpha / 8;
        float cos_theta = glm::dot(w, h);
        return cos_theta / (cos_theta * (1 - k) + k);
        // float out       = 0;
        // float cos_theta = glm::dot(w, h);
        // float x         = cos_theta > 0 ? 1 : 0;
        // float tan_theta = glm::tan(glm::acos(cos_theta));
        // float alpha     = 1 / (roughness * tan_theta);
        // float t         = roughness * tan_theta;
        ////beckmann
        // float A = alpha < 1.6 ? (1 - 1.259 * alpha + 0.396 * alpha * alpha) / (3.535 * alpha + 2.181 * alpha * alpha) : 0;
        // ggx
        // float A = (std::sqrt(1 + 1 / (alpha * alpha)) - 1) / 2;
        // return x / (1 + A);
    }

    float smith(float const & roughness, glm::vec3 const & wi, glm::vec3 const & w0, glm::vec3 const & h) {
        return smithG1(roughness, wi, h) * smithG1(roughness, w0, h);
    }

    glm::vec3 FR(float const & roughness, float const & shininess, glm::vec3 const & kd, glm::vec3 const & ks, glm::vec3 const & wi, glm::vec3 const & w0, glm::vec3 const & N) {
        // calculate the contribution of diffuse model
        glm::vec3 diffuse = glm::dot(N, w0) > 0.0f ? kd / PI : glm::vec3(0.0f);
        if (glm::dot(ks, ks) < EPS1) return diffuse;
        glm::vec3 h        = glm::normalize(wi + w0);
        float     f        = schlick(shininess, wi, N);
        float     g        = smith(roughness, wi, w0, N);
        float     d        = D(roughness, h, N);
        float     FGD      = f * g * d;
        glm::vec3 specular = ks * FGD / std::max((4 * glm::dot(wi, N) * glm::dot(w0, N)), 0.001f);
        return (1.0f - f) * diffuse + specular;
    }

    glm::vec3 SamplePoint(VCX::Engine::Light const & light) {
        static std::random_device                    dev;
        static std::mt19937                          rng(dev());
        static std::uniform_real_distribution<float> rand(0.0f, 1.0f);
        float                                        random_u = rand(rng);
        float                                        random_v = rand(rng);
        return light.Position + 130 * random_u * glm::vec3(1, 0, 0) + 115 * random_v * glm::vec3(0, 0, 1) - glm::vec3(65, 0, 57.5);
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

    glm::vec3 PathTrace(const RayIntersector & intersector, Ray ray, bool const & enableImprotanceSampling) {
        glm::vec3 color(0.0f);
        RayHit    rayHit = intersector.IntersectRay(ray);
        if (! rayHit.IntersectState) return color;
        const glm::vec3    rayDir    = glm::normalize(-ray.Direction);
        const glm::vec3    pos       = rayHit.IntersectPosition;
        const glm::vec3    n         = rayHit.IntersectNormal;
        const glm::vec3    kd        = rayHit.IntersectAlbedo;
        const glm::vec3    ks        = rayHit.IntersectMetaSpec;
        const float        alpha     = rayHit.IntersectAlbedo.w;
        const float        shininess = rayHit.IntersectMetaSpec.w * 256;
        VCX::Engine::Light light     = intersector.InternalScene->Lights[0];
        glm::vec3          lightNorm(0, -1, 0); // only for cornell box
        if (isLight(pos, light))
            return light.Intensity;
        float roughness = 0.2f;
        //对光源采样
        float     lightPdf    = 1.0f / (130 * 115);
        glm::vec3 lightPos    = SamplePoint(light);
        glm::vec3 l           = lightPos - pos;
        glm::vec3 lightDir    = glm::normalize(l);
        float     attenuation = 1.0f / glm::dot(l, l);
        auto      hit         = intersector.IntersectRay(Ray(lightPos, -l));
        float     dis         = glm::distance(pos, hit.IntersectPosition);
        while (hit.IntersectState && dis > EPS1) {
            if (hit.IntersectAlbedo.w >= 0.2f) break;
            hit = intersector.IntersectRay(Ray(hit.IntersectPosition, -l));
            dis = glm::distance(pos, hit.IntersectPosition);
        }
        if (! hit.IntersectState || dis <= EPS1) {
            glm::vec3 fr = FR(roughness, shininess, kd, ks, rayDir, lightDir, n);
            color += light.Intensity * fr * glm::dot(lightDir, n) * std::max(0.0f, glm::dot(-lightDir, lightNorm)) / lightPdf * attenuation;
        }
        //对其他方向
        static std::random_device                    dev;
        static std::mt19937                          rng(dev());
        static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        if (dist(rng) > RussianRoulette) return color;
        float     pdf     = 0.0f;
        glm::vec3 nextDir = glm::normalize(SampleLight(roughness, ks, rayDir, n, pdf, enableImprotanceSampling));
        Ray       nextRay(pos, nextDir);
        RayHit    nextRayHit = intersector.IntersectRay(nextRay);
        glm::vec3 nextPos    = nextRayHit.IntersectPosition;
        if (nextRayHit.IntersectState && ! isLight(nextPos, light)) {
            glm::vec3 fr = FR(roughness, shininess, kd, ks, rayDir, nextDir, n);
            color += PathTrace(intersector, nextRay, enableImprotanceSampling) * fr * glm::dot(nextDir, n) / pdf / RussianRoulette;
        }
        return color;
    }
} // namespace VCX::Labs::Rendering