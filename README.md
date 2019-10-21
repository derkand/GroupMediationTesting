# GroupMediationTesting
Authors: Andriy Derkach, Josh Sampson

Maintainer: Andriy Derkach (andriy.derkach@nih.gov)

Description: We consider the scenario where there is an exposure, multiple biologically-defined sets of biomarkers, and an outcome. We propose a new two-step procedure that tests if any of the sets of biomarkers mediate the exposure/outcome relationship, while maintaining a prespecified Family-Wise Error Rate (FWER). The first step of the proposed procedure is a screening step that removes all groups that are unlikely to be strongly associated with both the exposure and the outcome. The second step adapts recent advances in post-selection inference to test if there are true mediators in each of the remaining, candidate sets. Code for five tests are presented here.

Depndens: R (>= 3.2.1), psych, pscl, optimx, nleqslv, mvnfast, CompQuadForm
