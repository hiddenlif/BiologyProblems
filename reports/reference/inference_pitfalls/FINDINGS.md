# Executed findings: inference pitfalls

All data are synthetic and seeded. No patient, employer, proprietary assay or product data are included.

## Observed results

Target frequency increases from 20% to 35%, while the estimated target count falls from 200,000 to 140,000. Treating 96 technical pairs as independent produces p = 0.00104; the donor-level exact sign-flip p is 0.52344. In the fixed null simulation, random-row balanced accuracy is 1.000 versus 0.43125 with donor holdouts.

## Scope of evidence

The paired t interval assumes independent donors and an adequate model for donor differences. The sign-flip test assumes symmetric/exchangeable signs under the null. The deliberately unusual synthetic donor effects illustrate sensitivity to the unit of analysis, not a definitive population inference. The leakage example is one fixed negative control, not an estimate of typical inflation. Absolute counts assume representative flow sampling and an independently valid total viable-cell count.

The accompanying validation receipt identifies the analysis and dependency hashes and reports pytest and notebook execution. Re-run after changing code, sources or dependencies.
