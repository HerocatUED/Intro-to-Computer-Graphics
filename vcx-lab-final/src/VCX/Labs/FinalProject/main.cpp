#include "Assets/bundled.h"
#include "Labs/FinalProject/App.h"

int main() {
    using namespace VCX;
    return Engine::RunApp<Labs::Rendering::App>(Engine::AppContextOptions {
        .Title      = "VCX Labs Final: PathTracing",
        .WindowSize = { 1024, 768 },
        .FontSize   = 16,

        .IconFileNames = Assets::DefaultIcons,
        .FontFileNames = Assets::DefaultFonts,
    });
}
