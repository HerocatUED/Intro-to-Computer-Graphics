#pragma once

#include <vector>

#include "Engine/app.h"
#include "Labs/FinalProject/CasePathTracing.h"
#include "Labs/Common/UI.h"

namespace VCX::Labs::Rendering {
    class App : public Engine::IApp {
    private:
        Common::UI         _ui;

        CasePathTracing     _casePathTracing;

        std::size_t        _caseId = 0;

        std::vector<std::reference_wrapper<Common::ICase>> _cases = {         
            _casePathTracing,
        };

    public:
        App();

        void OnFrame() override;
    };
}
