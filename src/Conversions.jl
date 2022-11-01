"""
Ingestion rate (mg dwt/day) as a function of carapace length (cm) according to Urabe & Watanabe; originally determined for D. galeata.
"""
Urabe25(L) = 2*1e-3*12^1.75

"""
Length-dry weight relationship (cm vs mg dwt) according to Gelller & Mueller (1985).
"""
gellermass(L)=exp((1.5674+0.0287 * (60))+log(L) * (3.3611+0.0111 * 60)) / 10