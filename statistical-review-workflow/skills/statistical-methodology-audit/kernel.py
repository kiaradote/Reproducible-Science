
def design_gate(observational_unit, independent_sampling_unit, replication_unit, assumptions, source):
    """Establish the design before any methodological criticism. All arguments must be non-empty strings
    (assumptions may be a list). Returns the gate record; refuses empty values."""
    gate = {"observational_unit": observational_unit, "independent_sampling_unit": independent_sampling_unit,
            "replication_unit": replication_unit, "assumptions": assumptions, "source": source}
    missing = [k for k, v in gate.items() if not v]
    if missing:
        raise ValueError("design_gate incomplete: " + ", ".join(missing) + ". Do not criticise methods yet.")
    return gate


def concern_record(concern_id, claim_ids, error_type, statement, source_evidence, reasoning,
                   verification_procedure, falsification_attempt, result, confidence, impact,
                   evidence_tag, alternative_explanations=""):
    """One concern in the concern ledger.
    error_type: arithmetic / methodology / interpretation / reporting / generalization (consequences differ).
    result: demonstrated / supported concern / withdrawn / not assessable.
    'demonstrated' requires evidence_tag Verified, a performed verification_procedure and a falsification attempt.
    Every concern needs a falsification_attempt, including those left as 'not assessable'."""
    types = ("arithmetic", "methodology", "interpretation", "reporting", "generalization")
    results = ("demonstrated", "supported concern", "withdrawn", "not assessable")
    levels = ("low", "medium", "high")
    tags = ("Verified", "Inferred from text", "Not assessable")
    if error_type not in types:
        raise ValueError("error_type must be one of " + str(types))
    if result not in results:
        raise ValueError("result must be one of " + str(results))
    if confidence not in levels or impact not in levels:
        raise ValueError("confidence and impact must be low/medium/high")
    if evidence_tag not in tags:
        raise ValueError("evidence_tag must be one of " + str(tags))
    for name, val in (("statement", statement), ("source_evidence", source_evidence), ("reasoning", reasoning),
                      ("falsification_attempt", falsification_attempt)):
        if not val:
            raise ValueError(name + " is required")
    if result == "demonstrated" and (evidence_tag != "Verified" or not verification_procedure):
        raise ValueError("a demonstrated error needs evidence_tag=Verified and a performed verification_procedure")
    score = {"low": 1, "medium": 2, "high": 3}
    return {"concern_id": concern_id, "claim_ids": claim_ids, "error_type": error_type, "statement": statement,
            "source_evidence": source_evidence, "reasoning": reasoning,
            "verification_procedure": verification_procedure, "falsification_attempt": falsification_attempt,
            "alternative_explanations": alternative_explanations, "result": result, "confidence": confidence,
            "impact": impact, "evidence_tag": evidence_tag,
            "priority": score[confidence] * score[impact]}
