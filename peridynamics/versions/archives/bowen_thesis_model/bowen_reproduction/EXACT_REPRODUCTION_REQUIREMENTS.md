# Requirements to reproduce Bowen's original result

The objective is now **the original Chapter 3 computation**, rather than a physically corrected alternative. Do not infer that every numerical correction should be applied to the author's executable: reproducing a published program and improving its mechanics are separate tasks. Existing baseline/corrected outputs remain preserved.

## What is established

The thesis specifies geometry, many bulk parameters, mesh spacing and timestep, plus rounded aggregate outputs. It does not fully specify executable initialization and discretization. Supplied `Console2.f90` has dx19mm, dt3us, E30.5GPa, concrete density2500, 26mm bars and other loading; it is not verified as the Chapter3 solver. Assignment and line evidence is in `outputs/original_source_identity.json`. The alternate `Console322.f90` is a200×100mm two-dimensional model, also not Chapter3.

A further fidelity issue: the supplied Fortran's concrete force evaluation, like its interface evaluation, does not visibly consult previous failure flags. Therefore the current irreversible concrete/interface engine is not literal equivalence to that program. The omitted author Chapter3 context must resolve this before changing rupture history. Physical conservation of an alternative does not resolve author-code equivalence.

## Minimum author package

Prefer the exact original **source plus all input/initialization files and original output data** for §3.1 (38.4kg,10.8m/s; Table3-1 and Figure3-3). Include:

1. Initialization/assignment of `bcc,s01,sc1,vol` and steel/integration volumes; reinforcement node/connectivity arrays, virtual sample points, shape weights and family truncation.
2. Concrete force law and complete surrounding rate-code context: meaning/sign of `s1`, base threshold resets, AppendixA vs Eq2-13 empirical exponent, whether modulus is multiplied by DIF_E, prior `fail`/`failfail` handling and unloading.
3. Surface-energy correction implementation, quadrature and particle coordinate generation; concrete self-contact algorithm and neighbour search.
4. Hammer initial location, integration and recoil, rotation/guide constraints, contact ordering and force aggregation; support masks/widths, lift-off and anchors.
5. Rod constitutive/yield/rupture rules and full bent-end bar geometry.
6. Exact monitored displacement node/coordinate and deflection definition; original force/deflection time histories, final step, damage snapshots and output cadence.
7. Compiler/version/options, precision, thread count, any seeds, original executable if available, and observed author run outputs.

A Chinese request suitable for sending to the author (no message has been sent):

> 为复现硕士论文第三章3.1节的普通混凝土RC梁冲击计算（38.4kg、10.8m/s，表3-1的359.8kN和22.70mm），希望获取生成该结果的完整程序、初始化/输入文件以及原始冲击力、跨中位移和损伤输出。尤其需要钢筋—混凝土耦合参数bcc、s01、sc1、vol及控制体积、钢筋与耦合点坐标/连接关系、表面修正、接触与支座处理、应变率及断键历史处理、跨中测点定义和计算结束时间。请同时提供编译器版本与选项。现有Console2.f90使用19mm网格与不同构件/加载，尚无法确认为第三章原程序。

## Reproduction sequence after receipt

1. Archive and hash inputs; prove benchmark identity from assignments and generated geometry, without overwriting current work.
2. Run original program unmodified in a compatible runtime; compare raw outputs to author's originals. Diagnose platform differences before port changes.
3. Compare Python against Fortran at elementary interactions, complete one-step force/state, first contact, pre-peak, peak and unloading. Carry literal original behaviors explicitly; retain physically corrected variants separately.
4. Reproduce rounded Table3-1 values and original curves/diagnostic snapshots, with the author's exact measurement and output cadence. Preserve both every-step peaks and original reporting conventions.
5. Only claim original-results reproduction when the full run passes these comparisons. Mesh refinement and scientific validation are additional evidence, not replacements for matching the author's actual run.

The published values359.8kN and22.70mm are rounded. Without raw outputs and identical arithmetic/platform, **bitwise identity cannot be tested**. Initial acceptance is agreement at the precision printed in the thesis plus original history/crack morphology; no arbitrary material-factor tuning is permitted.

## Current blocker

The original Chapter3 executable definition is not available in inspected files. The unidentified interface values and source-history ambiguities cannot be derived uniquely from two peaks and a crack image. More runs of the provisional reconstruction cannot establish exact original-results reproduction. The user confirmed that only the supplied files are available. Do not ask again for the original author package. Proceed with source-supported reconstruction, preserve explicit unresolved choices, and do not promise exact author-run equivalence.
