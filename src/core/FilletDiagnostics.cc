#include "core/FilletDiagnostics.h"

#include <cassert>
#include <utility>

namespace {

thread_local unsigned int evaluationDepth = 0;
thread_local std::optional<std::string> pendingError;

}  // namespace

namespace FilletDiagnostics {

EvaluationScope::EvaluationScope()
{
  if (evaluationDepth++ == 0) pendingError.reset();
}

EvaluationScope::~EvaluationScope()
{
  assert(evaluationDepth > 0);
  evaluationDepth--;
}

void setError(std::string message)
{
  if (!pendingError) pendingError = std::move(message);
}

bool hasError()
{
  return pendingError.has_value();
}

std::optional<std::string> takeError()
{
  if (evaluationDepth > 0) return pendingError;
  auto result = std::move(pendingError);
  pendingError.reset();
  return result;
}

}  // namespace FilletDiagnostics
