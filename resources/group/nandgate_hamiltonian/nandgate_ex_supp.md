---
title: Encoding Logic Gates in Ising Hamiltonian
author: David Grisham
header-includes:
    \usepackage{braket}
    \newcommand{\su}{\ensuremath{\uparrow}}
    \newcommand{\sd}{\ensuremath{\downarrow}}
    \newcommand{\trans}{\ensuremath{\rightarrow}}
    \def\tightlist{}
...
\vspace{-0.8cm}

Determining Field ($h_i$) and Interaction ($J_{ij}$)Terms
=========================================================

General Outline of the Steps
----------------------------

- Start with any spin configuration, then start flipping individual qubits and
  think about how the energy of the system should *change* with those flips.
- When doing this, keep in mind the following:
    1. The lowest energy states should correspond to the states of the
       logic-gate of interest. So, when switching from a state that encodes one
       input/output set for our gate to one that does not, the energy should
       increase. We can think of this in terms of work/energy from classical
       mechanics, with this case being analogous to nonzero work being done on
       our system, since the total energy of the system changes.
    2. The states that correspond to our logic gate all have the same energy.
       This puts a restriction on how certain values (terms in the total energy
       equation) have to change with bit flips (i.e. do a bit flip, see which
       interaction terms and which field terms determine the change in energy of
       the system. If assumptions have been made about most of them, then the
       fact that the total energy should be the same can dictate the value of
       the final field/interaction terms we need). Following the work/energy
       analogy again, this would be the case where there is no work done in our
       system yet the relative values of internal energies change.


Concrete Example
================

As an example, let's look at how certain bit flips change the energy of our
system in the NAND gate example that Prof. Coffey provided, where A and B are
the inputs to our NAND gate and C is the output. Recall that spin UP (denoted
\su) corresponds to logical True, and spin DOWN (\sd) corresponds to logical
False. Assume the state of our system is written as $\ket{s_a s_b s_c}$, where
$s_a$, $s_b$, and $s_c$ correspond to the spins of inputs A, B, and C,
respectively. Also recall that the Ising Hamiltonian takes the following form:

$$H_f = \sum_j h_j s_j + \sum_{i < j} J_{ij} s_i s_j$$

Case 1: Non-NAND to NAND Transition
-----------------------------------

Consider how the system's energy changes in the following transition:

$$\ket{\sd\sd\sd} \trans \ket{\sd\sd\su}$$

Notice first that the starting state, $\ket{\sd\sd\sd}$, does not correspond to
one of the NAND gate configurations. Also notice that they state we end up in,
$\ket{\sd\sd\su}$, corresponds to A and B being False and C being True, which
agrees with the NAND gate logic. Because of this, we know that the energy of our
system should be lower in the second state, since the states that agree with the
NAND gate logic should be the only ground-states of our system. We also need to
pay attention to the spin that flipped in this transition, which was $s_c$. This
means that, when $s_c$ changed, the energy of our system decreased. There are
two key possible reasons for this:

1. $s_c$ flipped, which also means that its alignment with the $\vec{B}$-field
   term $h_c$ flipped as well. We know the energy decreased, so it's possible
   that it is now anti-aligned with the component of the $\vec{B}$-field that
   corresponds to this qubit. However, it is also possible that this term
   increases the total energy, since we currently do not have an alignment
   assigned to the $\vec{B}$-field and there are other terms we must consider
   that could have decreased the total energy as well. We can only say that
   *magnitude* of the change caused by this term is $2 \times h_c$.
2. Before the transition, $s_c$ was aligned with $s_a$ and $s_b$; after the
   transition, it was not. So the interaction constants $J_{ac}$ and $J_{bc}$
   are of interest here. We know that $J_{ac}=J_{bc}$, because A and B are
   equivalent as far as C is concerned. Let's assume that the interaction terms
   $J_{ac}$ and $J_{bc}$ are positive. Then, looking at the interaction terms in
   the Ising Hamiltonian, we see that if two spins $s_i$ and $s_j$ are aligned
   and their interaction term $J_{ij}$ is positive, then the energy of the
   system includes a **positive** $J_{ij}$ term. If they are anti-aligned, the
   total energy includes a **negative** $J_{ij}$ term. Given that C went from
   being aligned with A and B to being anti-aligned with them, the total energy
   of our system must have *decreased* by $2 \times (J_{ac} + J_{bc})$.

Now we can write an equation that relates the non-zero changes in energy within
our system to the overall change (which, at this point, we only know must be
less than 0):

$$\pm 2 \times h_c - 2 \times (J_{ac} + J_{bc}) < 0$$

Notice that we still cannot say anything certain about the sign of the field
term, despite the assumptions made about the signs of the interaction terms.

We made one assumption above, namely that the interaction terms $J_{ac}$ and
$J_{bc}$ were positive. This was not necesary, as there are two possible terms
that could decrease the energy of our system. However, making assumptions such
as this narrows down the form that our problem will take; we could have said
that those terms were negative instead, which would change the arguments we
could make in subsequent parts of the analysis.

Case 2: NAND to NAND Transition
-------------------------------

Now let's consider a second transition, one that goes between two states that
agree with the NAND gate logic (and, thus, each correspond to a ground state of
our system):

$$\ket{\sd\sd\su} \trans \ket{\sd\su\su}$$

This time, we flipped $s_b$. We know from the previous part that the first of
these states agrees with the NAND logic, and we can verify that the second
state, $\ket{\sd\su\su}$, does as well (A is False, B is True, so C should be
True as well). As mentioned previously, that fact that both states satisfy our
gate logic means that they should both correspond to the same ground state
energy. We know that our field and interaction constants are not simply 0, so
there must have been a tradeoff between certain energy values internal to our
system. Specifically, for the transition in question:

1. $s_b$ flipped, which means that its alignment with the relevant
   $\vec{B}$-field term, $h_b$, must have flipped as well. As in the first case,
   we cannot tell without further analysis whether this term increased or
   decreased during the transition. The magnitude of the change caused by this
   term is $2 \times h_b$.
2. Before the transition $s_b$ was aligned with $s_a$ only, and after the
   transition it was aligned with only $s_c$. Our previous assumption that
   $J_{bc}$ is positive tells us that the corresponding term in the total energy
   equation *increased* by $2 \times J_{bc}$. However, we still do not have
   enough information to tell whether $J_{ab}$ is positive or negative. However,
   since we already assumed $J_{ac}$ and $J_{bc}$ were positive, we can simply
   argue for consistency and say that $J_{ab}$ is also positive. This means that
   corresponding term in our total energy equation *decreased* by $2 \times
   J_{ab}$, since $s_a$ is spin DOWN and $s_b$ went from spin DOWN to spin UP
   (which caused the sign flip).

Given all of this information, we can now characterize the relationship between
the energy changes within our system as follows:

$$\pm 2 \times h_b + 2 \times (J_{bc} - J_{ab}) = 0$$


Concluding Thoughts
===================

The process presented here is a rough outline of how we might go about assigning
signs and values to the field and interaction terms for a Hamiltonian that
represents a logic gate. There are steps that could be taken, such as looking at
additional transitions and the corresponding changes in energy, to further
simplify the analysis above. Once an adequate number of constraints have been
placed on the variables, one could start assigning signs (as was done above) and
concrete values to the constants, since we really just care about the relative
signs and magnitudes of the $h_i$ and $J_{ij}$ values. Thus, there are multiple
sets of values that could be assigned to these constants, as long as they obey
all of the necessary constraints.

