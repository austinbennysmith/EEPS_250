## Coupling Daisy World to Lotka Volterra!
This is where I would write the report. Above would be my code. I would hyperlink to code throughout the report to make it easy to see what code produced which graphs. For example, [here is a link to the Lotka Volterra code I've uploaded already](https://github.com/austinbennysmith/EEPS_250/blob/main/Lotka_alone.m)

One unrealistic aspect of LV models outlined here: https://vlab.amrita.edu/?sub=3&brch=67&sim=185&cnt=1
"Each predator eats a constant proportion of the prey population per year; In other words, doubling the prey population will double the number eaten per predator, regardless of how big the prey population is" - no satiation

Show that 3-level LV model isn't working due to some kind of iteration issue.

A few sources:
http://math.bd.psu.edu/faculty/jprevite/mathmag243-255.pdf
https://sites.math.washington.edu/~morrow/336_16/2016papers/lalith.pdf

https://link.springer.com/article/10.1007%2Fs12043-020-1942-9


In the 3 level program: the reason the cycle slows down over time is that the herbivore death rate is ONLY tied to the number of carnivores, so if the number of carnivores crashes low enough the recovery is extremely slow

G=0.9 is more helpful than G=5 in the level 3 program because it allows one to see the carnivore fluctuations. Or at least lower G. See if lower G works for level 2 as well for consistency. It seems that in level 3 in the range of luminosities creating a Lotka-Volterra pattern, the time scale can be increased to *10^4 to get a pattern that is qualitatively similar to that observed in level 2 for the Lotka-Volterra-esque luminosities. However there is a time lag and a slowdown, and showing time at *10^1 time scale may be helpful to demonstrate this.
