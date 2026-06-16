# OpenSMA: Clinical PBPK and Systems Biology Mathematical Modeling Guide

This document presents the complete mathematical derivations and differential equations (ODEs) for the multi-compartment physiologically-based pharmacokinetics (PBPK) model, pharmacodynamics (PD), and Markov cohort simulations used in the OpenSMA platform.

---

## 1. Multi-Compartment PBPK Differential Equations (ODEs)

The platform simulates drug absorption, distribution, metabolism, and excretion (ADME) using a 5-compartment continuous mass transfer ODE system.

### 1.1. Compartment Definitions and Volume Scaling
The compartments are defined below:
1.  **Gut ($V_{\text{gut}}$):** Site of oral administration (for small molecules).
2.  **Plasma ($V_{\text{plasma}}$):** Systemic circulation.
3.  **Tissue ($V_{\text{tissue}}$):** Peripheral tissues (liver, kidneys, muscles).
4.  **CSF ($V_{\text{csf}}$):** Cerebrospinal fluid. Site of intrathecal bolus injection (for ASOs).
5.  **Brain Parenchyma ($V_{\text{brain}}$):** Intracellular site of action.

Compartmental volumes ($V$, mL) are scaled dynamically based on patient age ($t_{\text{age}}$, months):

$$V_i(t_{\text{age}}) = V_{i\text{, adult}} \times \left( \frac{\text{Weight}(t_{\text{age}})}{\text{Weight}_{\text{adult}}} \right)$$

### 1.2. The Mass Transfer ODE System
The concentration kinetics ($\mu\text{g/mL}$) for each compartment are governed by the following coupled differential equations:

#### 1. Gut Compartment (Oral Intake)
$$\frac{dC_{\text{gut}}}{dt} = -k_a \cdot C_{\text{gut}}$$

Where $k_a$ ($h^{-1}$) is the absorption rate constant (defined as $0.0$ for intrathecally injected ASOs).

#### 2. Plasma Compartment
$$\frac{dC_{\text{plasma}}}{dt} = \frac{k_a \cdot V_{\text{gut}} \cdot C_{\text{gut}}}{V_{\text{plasma}}} + k_{\text{diff\_tissue}} \cdot (C_{\text{tissue}} - C_{\text{plasma}}) + k_{\text{diff\_bbb}} \cdot \frac{V_{\text{brain}}}{V_{\text{plasma}}} \cdot (C_{\text{brain}} - C_{\text{plasma}}) - k_{\text{elim}} \cdot C_{\text{plasma}}$$

Where:
-   $k_{\text{diff\_tissue}}$: Mass exchange coefficient with peripheral tissue.
-   $k_{\text{diff\_bbb}}$: Active/passive permeability coefficient across the Blood-Brain Barrier (BBB).
-   $k_{\text{elim}}$: Renally/hepatically driven systemic clearance rate.

#### 3. Tissue Compartment
$$\frac{dC_{\text{tissue}}}{dt} = k_{\text{diff\_tissue}} \cdot \frac{V_{\text{plasma}}}{V_{\text{tissue}}} \cdot (C_{\text{plasma}} - C_{\text{tissue}})$$

#### 4. Cerebrospinal Fluid (CSF) Compartment (Intrathecal Injection Site)
ASOs are administered directly into this compartment as a bolus dose. Mass transport to the brain parenchyma and clearance via CSF bulk flow:
$$\frac{dC_{\text{csf}}}{dt} = k_{\text{diff\_csf}} \cdot \frac{V_{\text{brain}}}{V_{\text{csf}}} \cdot (C_{\text{brain}} - C_{\text{csf}}) - k_{\text{clear\_csf}} \cdot C_{\text{csf}}$$

#### 5. Brain Parenchyma Compartment (Target Site of Action)
$$\frac{dC_{\text{brain}}}{dt} = k_{\text{diff\_bbb}} \cdot (C_{\text{plasma}} - C_{\text{brain}}) + k_{\text{diff\_csf}} \cdot (C_{\text{csf}} - C_{\text{brain}}) - k_{\text{metab\_brain}} \cdot C_{\text{brain}}$$

Where $k_{\text{metab\_brain}}$ represents local metabolic degradation of the therapeutic agent inside brain cells.

---

## 2. Pharmacodynamics (PD) and Splicing Kinetics

Local therapeutic concentration ($C_{\text{local}}(t) \equiv C_{\text{brain}}(t)$) modulates the pre-mRNA splicing machinery, increasing the synthesis rate of functional SMN protein.

### 2.1. Splicing Modulation via Hill Equation
The fold-increase in SMN protein expression ($Boost_{\text{SMN}}$) is modeled using a cooperative Hill equation:

$$\text{Boost}_{\text{SMN}}(t) = \frac{E_{\max} \cdot C_{\text{local}}(t)^n}{EC_{50}^n + C_{\text{local}}(t)^n}$$

Where:
-   $E_{\max}$: Maximum asymptotic splicing induction fold-increase.
-   $EC_{50}$: Target concentration required to achieve 50% of $E_{\max}$ (affinity proxy).
-   $n$: Hill cooperativity coefficient ($n = 1.2$ in our calibrated models).

Total functional SMN protein concentration ($smn(t)$) relative to baseline:

$$smn(t) = smn_{\text{baseline}} \cdot (1.0 + \text{Boost}_{\text{SMN}}(t))$$

### 2.2. Motor Neuron Survival Dynamics
In untreated Type 1 SMA, the motor neuron population ($MN$) undergoes rapid apoptotic degeneration. Restored SMN protein ($smn(t)$) slows down this decay:

$$\text{protection}(t) = \min\left(1.0, \max\left(0.0, \frac{smn(t)}{100.0}\right)\right)$$

$$\frac{d[MN]}{dt} = -r_{\text{decay}} \cdot \left(1.0 - \text{protection}(t) \cdot \eta \right) \cdot [MN]$$

Where:
-   $r_{\text{decay}}$: Basal motor neuron degeneration rate ($0.05 \text{ month}^{-1}$) in untreated Type 1 SMA.
-   $\eta$: Splicing protection efficiency ($0.85$ to $0.95$).
-   $[MN]$: Percentage of viable motor neurons relative to birth healthy baseline.

---

## 3. Markov Modeling and Monte Carlo Cohort Simulation

To project clinical outcomes for a virtual cohort of 10,000 patients, we map the motor neuron survival fraction ($[MN](t)$) to clinical motor milestone transitions using a Markov chain.

### 3.1. Milestone Probability Functions
The probability of gaining key motor milestones (sitting, standing, walking) is modeled as a non-linear sigmoid function of motor neuron pool size:

$$P(\text{Sit} \mid [MN]) = \frac{1.0}{1.0 + e^{-0.1 \cdot ([MN] - 40.0)}}$$
$$P(\text{Stand} \mid [MN]) = \frac{1.0}{1.0 + e^{-0.1 \cdot ([MN] - 55.0)}}$$
$$P(\text{Walk} \mid [MN]) = \frac{1.0}{1.0 + e^{-0.1 \cdot ([MN] - 70.0)}}$$

### 3.2. Markov Transition Matrix
At each time step $\Delta t = 1$ month, the state of each patient transitions according to a transition probability matrix $M(t)$:

$$M(t) = \begin{pmatrix} 
P_{D\to D} & P_{D\to MN\_loss} & 0 & 0 \\
0 & P_{MN\_loss\to MN\_loss} & P_{MN\_loss\to Sit} & 0 \\
0 & 0 & P_{Sit\to Sit} & P_{Sit\to Walk} \\
0 & 0 & 0 & 1.0 
\end{pmatrix}$$

Where:
*   **Deceased (D):** Terminal state representing severe motor neuron depletion ($[MN] < 10\%$).
*   **State transitions** are constrained by the cumulative probability functions computed from $[MN](t)$ at each step.
*   **Survival Probability:** Estimated via Kaplan-Meier estimator over the simulated cohort:
    $$\hat{S}(t) = \prod_{t_i \le t} \left( 1 - \frac{d_i}{n_i} \right)$$
    Where $d_i$ is the number of deaths at time $t_i$, and $n_i$ is the number of patients surviving up to $t_i$.
