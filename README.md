# Exposing the Variety of Nash-Equilibria in Oligopolistic Electricity Market

Abstract
--------
Liberalized electricity markets promise a cost-efficient operation and expansion of power systems but may
as well introduce opportunities for strategic gaming for price-making agents. Given the rapid transition
of today’s energy systems, unconventional generation and consumption patterns are emerging, presenting
new challenges for regulators and policymakers to prevent such behaviors. The strategic offering of various
price-making agents in oligopolistic electricity markets resembles a multi-leader-common-follower game. The
decision problem of each agent can be modeled as a bi-level optimization problem, consisting of the strategic
agent’s decision problem in the upper-level, and the market clearing problem in the lower-level. When
modeling a multi-leader game, i.e., a set of bi-level optimization problems, the resulting equilibrium problem
with equilibrium constraints poses several challenges. Real-life applicability or policy-oriented studies are
challenged by the potential multiplicity of equilibria and the difficulty of exhaustively exploring this range of
equilbria. In this paper, the range of equilibria is explored by using a novel simultaneous solution method.
The proposed solution technique relies on applying Scholtes’ regularization before concatenating the strategic
actor’s decision problems’ optimality conditions. Hence, the attained solutions are stationary points with
high confidence. In a stylized example, different strategic agents, including an energy storage system,
are modeled to capture the asymmetric opportunities they may face when exercising market power. Our
analysis reveals that these models’ outcomes may span a broad range, impacting the derived economic metrics
significantly
This code is published as a companion material of the above named research, which is currently submitted to the European Journal of Operational Research and can be found on: https://www.mech.kuleuven.be/en/tme/research/energy-systems-integration-modeling/pdf-publications/wp-esim2020-02

Section 4.3. Convergence analysis in case of a single strategic agent
--------
The ``MILP_quantity.jl``, ``MPPDC_quantity.jl``, ``NLP_dual_MPPDC_quantity.jl`` files were used in this Section. All models can represent a single strategic actor (ESS,CG,RG) participating on the LL market, as explained in Section 4 of the working paper. An example usage of all three models is shown by ``test_script_singleagent.jl``.

Section 4.4. Numerical results of the stylized EPEC setting and 4.5. Analyzing the impact of changes in the strategic agents’ capacity
--------

``EPEC_quantitity.jl`` and ``diagonalization.jl`` files were used. An example usage of the EPEC model by changing the central planner's objective is shown by ``test_script_multiagent.jl``

License
--------
This work is licensed under a [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/)

Keywords
--------
OR in energy, Multi-leader-common-follower games, Electricity Market, Equilibrium Problem
with Equilibrium Constraints, Oligopoly
