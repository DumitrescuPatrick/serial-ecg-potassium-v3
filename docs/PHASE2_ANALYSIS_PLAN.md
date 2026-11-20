# Phase2 exploratory analysis plan —20November2026 UTC

Written before V3 amplitude-ratio scores or effect estimates are calculated in this run. A separate V3 Phase2 plan was not located in project or bounded local research-document searches; V1/V2 and earlier V3 feasibility decisions have already been examined. This plan is NOT preregistered or independent confirmation. Preserve every existing V3 file and the1672-record pre-completion manifest before work.

## Preconditions

Independently verify every saved pre-completion record against the full source index and recorded transfer hashes: subject/study/path/time,readable official-format headers,data length,units/calibration,digital checksums. Recheck all2037selected ECGs with the unchanged V3 engineering/extraction conditions,header time and local-copy checks. Preserve fixed selected pair IDs,cohort rules,age/encounter labels,±30min lab matches,2–48h pair interval and one pair/person/group. Verify actual extracted-RR shared counts are122/650hyper and200/650hypo; do not accept a stated number without checking. Stop effects if integrity or gate fails. The exact100/100gate and original≥3perarm shared-cell requirement remain unchanged.

## Endpoint and estimand

Use the original eight independent leads I,II,V1–V6,unchanged filter/QRS/post-QRS extraction and windows. Do not tune features,repair rules,exclusions,models,windows or seeds using V3 effects. L=log(post-QRS RMS/QRS RMS). Preserve historical numerical implementation: case scores from storedfloat32templates with1e-12RMS stabilizers; control scores from original double-precision feature RMS without the stabilizer. This small historical precision distinction is disclosed,not silently harmonized. Private worker instrumentation may retain original RMS values but must not change QC decisions or extracted RR.

Case signed change is abnormal minus normal regardless of chronological ordering; control signed change is later minus earlier. Primary outcome is absolute change. One pair/patient/group; cases may overlap hyper/hypo,controls disjoint. Primary effect is control-count-weighted difference in mean absolute change within the original interval bins[2,12),[12,24),[24,48]h × |deltaRR|≤100,(100,200],>200ms,requiring≥3patients in each arm/cell. Use EXACT V1 effect() function,seed20261119,2000patient resamples per arm/stratum,fixed observed control weights. Report estimate,percentile95%CI,eligible/common-support cases/controls,cells and missingness. Intervals condition on feature/cohort/strata and omit selection/training uncertainty. Shared controls/patients mean comparisons are not independent replications. No causal,diagnostic,clinical-utility or novel-phenotype claim.

## Fixed analysis matrix — report every row for both groups

1 Primary frozen30min cohort.
2 Both retained lab gaps≤15min; subset only,no rematching or pair replacement.
3 ±60min rule on frozen30min patients/matches: identical by construction. Explicitly mark not an independently rebuilt60min sensitivity; a newly linked60min cohort is OUT OF SCOPE because it would change the frozen design.
4 Both collections at/before their corresponding ECG; mixed direction omitted with counts.
5 Both collections at/after their corresponding ECG; simultaneous pairs can belong to both direction rows and must be disclosed.
6 Exclude any hemolysis-comment mention in either ECG’s lab,from BOTH arms; unadjudicated comments,not an assay result.
7 Exclude any index/machine,historical/header time disagreement from BOTH arms. V3 already excludes source disagreement; report identical results if no changes. Do not redefine clock verification.

Every row has eligible and common-stratum distinct-patient counts,effect/CI,or explicit insufficiency reason. No outcome-driven subgroup,covariate ormethod selection. Supplemental descriptive signed medians and chronology/age/sex/care/interval/timing/QC counts may be reported without inferential discovery claims. No P-values required; sensitivities are exploratory and not multiplicity-adjusted confirmation. Do not emphasize only positive rows.

## Execution and reporting

Freeze source/input/plan hashes,save an atomic completion receipt to prevent duplicate effect runs. Reuse verified outputs if the analysis already completed with identical hashes. No re-extraction/reselection after effect inspection. Numerical validation may independently recompute the prespecified statistic without tuning. Correct stale retrieval/gate text in a new superseding report; preserve old files as historical snapshots. Draft a full exploratory Cureus-format article only from verified outputs,aggregate figures/tables,limitations and declaration placeholders. Public ZIP must exclude IDs,dates,individual scores,paths,signals,templates,weights and review keys. Physical clock alignment,treatment timing,hemolysis adjudication,source release versions and independent clinical waveform interpretation remain unresolved unless actually verified.
