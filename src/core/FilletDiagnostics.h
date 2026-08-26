#pragma once

#include <optional>
#include <string>

namespace FilletDiagnostics {

class EvaluationScope
{
public:
  EvaluationScope();
  ~EvaluationScope();

  EvaluationScope(const EvaluationScope&) = delete;
  EvaluationScope& operator=(const EvaluationScope&) = delete;
};

void setError(std::string message);
bool hasError();
std::optional<std::string> takeError();

}  // namespace FilletDiagnostics
