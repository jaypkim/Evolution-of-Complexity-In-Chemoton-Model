extensions [cf csv ]
__includes ["collision.nls"
            "reaction.nls"
            "molecule-char.nls"
            "enzyme-kinetics.nls"
            "performance-rules.nls"
            "genome.nls"
            "data-collection.nls"
            "molecule-char-run.nls"
            ]

;;Global variables that affects all agents
globals [
  tick-delta                      ;; how much we advance the tick counter this time through
  max-tick-delta                  ;; the largest tick-delta is allowed to be
  init-avg-speed
  init-avg-energy
  init-speed
  viscosity
  avg-speed
  avg-energy
  total-energy
  init-substrate-conc
  performance
  mutation-prob
  performance-reproduce
  generation-reproduce
  performance-reproduce-record
  generations-reproduce-record
  generations
  implicit-prob-gone
  implicit-prob-here
  implicit-prob-gone-a-1-n
  implicit-prob-here-a-1-n
  implicit-prob-gone-a-2-n
  implicit-prob-here-a-2-n
  implicit-prob-gone-a-3-n
  implicit-prob-here-a-3-n
  implicit-prob-gone-a-4-n
  implicit-prob-here-a-4-n
  implicit-prob-gone-a-5-n
  implicit-prob-here-a-5-n
  implicit-prob-gone-a-1-s
  implicit-prob-here-a-1-s
  implicit-prob-gone-a-2-s
  implicit-prob-here-a-2-s
  implicit-prob-gone-a-3-s
  implicit-prob-here-a-3-s
  implicit-prob-gone-a-4-s
  implicit-prob-here-a-4-s
  implicit-prob-gone-a-5-s
  implicit-prob-here-a-5-s
  implicit-prob-gone-a-1-r
  implicit-prob-here-a-1-r
  implicit-prob-gone-a-2-r
  implicit-prob-here-a-2-r
  implicit-prob-gone-a-3-r
  implicit-prob-here-a-3-r
  implicit-prob-gone-a-4-r
  implicit-prob-here-a-4-r
  implicit-prob-gone-a-5-r
  implicit-prob-here-a-5-r
  implicit-prob-gone-mod-a-1-n
  implicit-prob-here-mod-a-1-n
  implicit-prob-gone-mod-a-2-n
  implicit-prob-here-mod-a-2-n
  implicit-prob-gone-mod-a-3-n
  implicit-prob-here-mod-a-3-n
  implicit-prob-gone-mod-a-4-n
  implicit-prob-here-mod-a-4-n
  implicit-prob-gone-sub-s
  implicit-prob-here-sub-s

  birth-mut-counter
  a-1-n-evo-prob
  a-2-n-evo-prob
  a-3-n-evo-prob
  a-4-n-evo-prob
  a-5-n-evo-prob
  a-1-s-evo-prob
  a-2-s-evo-prob
  a-3-s-evo-prob
  a-4-s-evo-prob
  a-5-s-evo-prob

  a-1-s-performance-1
  a-2-s-performance-1
  a-3-s-performance-1
  a-4-s-performance-1
  a-5-s-performance-1

  a-1-s-performance-2
  a-2-s-performance-2
  a-3-s-performance-2
  a-4-s-performance-2
  a-5-s-performance-2


  a-1-n-performance-1
  a-2-n-performance-1
  a-3-n-performance-1
  a-4-n-performance-1
  a-5-n-performance-1

  a-1-r-performance-1
  a-2-r-performance-1
  a-3-r-performance-1
  a-4-r-performance-1
  a-5-r-performance-1
]

breed [substrate substrates]
breed [a-1-n]
breed [a-1-s]
breed [a-2-n]
breed [a-2-s]
breed [a-3-n]
breed [a-3-s]
breed [a-4-n]
breed [a-4-s]
breed [a-5-n]
breed [a-5-s]
breed [a-1-r]
breed [a-2-r]
breed [a-3-r]
breed [a-4-r]
breed [a-5-r]
breed [mod-a-1-n]
breed [mod-a-2-n]
breed [mod-a-3-n]
breed [mod-a-4-n]
breed [genome]
breed [mod-a-1-s]
breed [mod-a-2-s]
breed [mod-a-3-s]
breed [mod-a-4-s]
breed [substrate-s]

turtles-own
[
  compound-name
  mass
  radius
  speed
  energy
  last-collision
  rxn-reactant
  rxn-type
  rxn-time-dif
  rxn-prob
  k_1
  k_-1
  performance-1
  performance-2
  rxn-tot-1
  rxn-tot-2
  mut-prob-run
  evo-tag
  r-tag
  next-gen-tag
  mutation-gamma-dist-alpha
  mutation-gamma-dist-lambda
]

to setup-globals
  ;; According to http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1225847/pdf/biophysj00084-0290.pdf
  ;; viscosity of cytoplasm is roughly 0.0011 kPa
  set viscosity 0.0011
  ;; Arbitrary value to allow tick progression to occur in a reasonable manner
  set max-tick-delta 0.1073
  set init-substrate-conc 100
  ;; if changed, change in the go procedure as well
  set performance 500
  set performance-reproduce 0
  set generation-reproduce 40
  set performance-reproduce-record 0
  set generations-reproduce-record 0
  set birth-mut-counter 0
  implicit-probability
  genome-evo
  genome-r
  resize-world -100 100 -100 100
end


to diffusion
    ;; 10 ^ 20 is an arbitrary number that allows the simulation to be seen,
    ;; as diffusion cannot move a measurable degree in meters
    set init-speed sqrt (2 * (1.3806495279 * 10 ^ -19 * temperature * 10 ^ 20) /
     (6 * pi * viscosity * radius))
end

;; Means of having appropriate speeds
to kinetic
  set energy (0.5 * mass * (speed ^ 2))
end

to setup
  clear-all
  reset-ticks
  setup-globals
  number-of-agents
  ask turtles [
    diffusion
    set rxn-reactant nobody
    set evo-tag []
    set r-tag []
  ]
  calculate-tick-delta
  substrate-char
  a-1-n-char
  a-2-n-char
  a-3-n-char
  a-4-n-char
  a-5-n-char
  mod-a-1-n-char
  mod-a-2-n-char
  mod-a-3-n-char
  mod-a-4-n-char
  mod-a-1-s-char-i
  mod-a-2-s-char-i
  mod-a-3-s-char-i
  mod-a-4-s-char-i
  substrate-s-char-i
  a-1-r-char
  a-2-r-char
  a-3-r-char
  a-4-r-char
  a-5-r-char
  a-1-s-char
  a-2-s-char
  a-3-s-char
  a-4-s-char
  a-5-s-char
  genome-char
  global-mutator-means
end

;;; Have to insert every rxn that is possible in the go function ;;;

to go
  ;; *Place end conditions here as a conditional*
  ask turtles
  [check-for-collision
  diffusion
  move
  if not (breed = genome) [set next-gen-tag 0]

  ;; Add rxn-functions here. For speed, write an if function
  ;; specifying the breed and rxn-reactant conditions, as shown below

  if breed = a-1-n and (rxn-reactant = "complex")
    [ask self [dissociate-a1-&-s 15 3 3]]
  if breed = a-1-n and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-1-n-&-sub-rxn]

  if breed = a-2-n and (rxn-reactant = "complex")
    [ask self [dissociate-a2-&-mod-a1 15 3 3]]
  if breed = a-2-n and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-2-n-&-mod-a-1-n-rxn]

  if breed = a-3-n and (rxn-reactant = "complex")
    [ask self [dissociate-a3-&-mod-a2 15 3 3]]
  if breed = a-3-n and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-3-n-&-mod-a-2-n-rxn]

  if breed = a-4-n and (rxn-reactant = "complex")
    [ask self [dissociate-a4-&-mod-a3 15 3 3]]
  if breed = a-4-n and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-4-n-&-mod-a-3-n-rxn]

  if breed = a-5-n and (rxn-reactant = "complex")
    [ask self [dissociate-a5-&-mod-a4 15 3 3]]
  if breed = a-5-n and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-5-n-&-mod-a-4-n-rxn]

;;;;;;;;;;;;;;
;; a-#-s
;;;;;;;;;;;;;;

  if breed = a-1-s and (rxn-reactant ="complex") and (rxn-type = 1)
  [ask self [dissociate-a1s-&-mod-a2 15 3 3]]
  if breed = a-1-s and (rxn-reactant ="complex") and (rxn-type = 2)
  [ask self [dissociate-a1s-&-substrate 15 3 3]]
  if breed = a-1-s and (rxn-reactant = nobody)
  [a-1-s-&-mod-a-2-n-rxn
    a-1-s-&-sub-rxn]

  if breed = a-2-s and (rxn-reactant ="complex") and (rxn-type = 1)
  [ask self [dissociate-a2s-&-mod-a3 15 3 3]]
  if breed = a-2-s and (rxn-reactant ="complex") and (rxn-type = 2)
  [ask self [dissociate-a2s-&-mod-a1 15 3 3]]
  if breed = a-2-s and (rxn-reactant = nobody)
  [a-2-s-&-mod-a-3-n-rxn
    a-2-s-&-mod-a1-rxn]

  if breed = a-3-s and (rxn-reactant ="complex") and (rxn-type = 1)
  [ask self [dissociate-a3s-&-mod-a4 15 3 3]]
  if breed = a-2-s and (rxn-reactant ="complex") and (rxn-type = 2)
  [ask self [dissociate-a3s-&-mod-a2 15 3 3]]
  if breed = a-2-s and (rxn-reactant = nobody)
  [a-3-s-&-mod-a-4-n-rxn
    a-3-s-&-mod-a2-rxn]

  if breed = a-4-s and (rxn-reactant ="complex") and (rxn-type = 1)
  [ask self [dissociate-a4s-&-sub 15 3 3]]
  if breed = a-4-s and (rxn-reactant ="complex") and (rxn-type = 2)
  [ask self [dissociate-a4s-&-mod-a-3 15 3 3]]
  if breed = a-2-s and (rxn-reactant = nobody)
  [a-4-s-&-sub-rxn
    a-4-s-&-mod-a3-rxn]

  if breed = a-5-s and (rxn-reactant = "complex") and (rxn-type = 1)
  [ask self [dissociate-a5s-&-mod-a-4 15 3 3]]
  if breed = a-5-s and (rxn-reactant = "complex") and (rxn-type = 2)
  [ask self [dissociate-a5s-&-mod-a-1 15 3 3]]
  if breed = a-5-s and (rxn-reactant = nobody)
  [a-5-s-&-mod-a4-rxn
    a-5-s-&-mod-a1-rxn]

;;;;;;;;;;
;; a-#-r
;;;;;;;;;;

;; Add rxn-functions here. For speed, write an if function
  ;; specifying the breed and rxn-reactant conditions, as shown below

  if breed = a-1-r and (rxn-reactant = "complex")
    [ask self [dissociate-a1r-&-s-s 15 3 3]]
  if breed = a-1-r and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-1-r-&-sub-s-rxn]

  if breed = a-2-r and (rxn-reactant = "complex")
    [ask self [dissociate-a2r-&-mod-a1s 15 3 3]]
  if breed = a-2-r and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-2-r-&-mod-a-1-s-rxn]

  if breed = a-3-r and (rxn-reactant = "complex")
    [ask self [dissociate-a3r-&-mod-a2s 15 3 3]]
  if breed = a-3-r and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-3-r-&-mod-a-2-s-rxn]

  if breed = a-4-r and (rxn-reactant = "complex")
    [ask self [dissociate-a4r-&-mod-a3s 15 3 3]]
  if breed = a-4-r and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-4-r-&-mod-a-3-s-rxn]

  if breed = a-5-r and (rxn-reactant = "complex")
    [ask self [dissociate-a5r-&-mod-a4s 15 3 3]]
  if breed = a-5-r and (rxn-reactant = nobody)
  ;; If an agent has more than 1 possible reaction, add it inside the below bracket
  [a-5-r-&-mod-a-4-s-rxn]

  if breed = genome
  ;; Every 3 ticks, genome produce is called
  [if ticks mod 3 = 0 [genome-produce]]

  ;; Mutation probability before hatching in genome. This command overrides the command
  ;; from the mutation command prior to a hatched molecule's data.

  ]


  ;calculate-tick-delta

  ;; Survival conditions
  if performance < -2000 or ticks > 2000 [
    stop
    ]
  ;; For the reproduction of the cell
  ;; induce mutation prob change after reproduction...?
  if (performance-reproduce - performance + performance-reproduce-record < 0)
  and
  (generation-reproduce - ticks + generations-reproduce-record < 0)
  [set generations (1 + generations)
    set performance-reproduce-record performance
    set generations-reproduce-record ticks
    set birth-mut-counter 1
    ;; if changed here, change above in setup as well
    ;; This performance change represents the new generation's ability
    ;; to handle the environment under specified RTPV
    set performance 500

    ifelse count a-1-n > 30[
    ask n-of 30 a-1-n
    [set next-gen-tag 1]]
    [ask a-1-n
      [set next-gen-tag 1]]

    ifelse count a-1-s > 30[
    ask n-of 30 a-2-n
    [set next-gen-tag 1]]

    [ask a-2-n
      [set next-gen-tag 1]]
    ifelse count a-3-n > 30[
    ask n-of 30 a-3-n
    [set next-gen-tag 1]]

    [ask a-3-n
      [set next-gen-tag 1]]
    ifelse count a-4-n > 30
    [
    ask n-of 30 a-4-n
    [set next-gen-tag 1]]

    [ask a-4-n
      [set next-gen-tag 1]]
    ;; a-5-n
    ifelse count a-5-n > 30[
    ask n-of 30 a-5-n
    [set next-gen-tag 1]]

    [ask a-5-n
      [set next-gen-tag 1]]
    ;; a-1-s
    ifelse count a-1-s > 30
    [ask n-of 30 a-1-s
      [set next-gen-tag 1]]

    [ask a-1-s
      [set next-gen-tag 1]]
    ifelse count a-2-s > 30
    [ask n-of 30 a-2-s
      [set next-gen-tag 1]]
    [ask a-2-s
      [set next-gen-tag 1]]
    ;; a-3-s
    ifelse count a-3-s > 30
    [ask n-of 30 a-3-s
      [set next-gen-tag 1]]
    [ask a-3-s
      [set next-gen-tag 1]]
    ifelse count a-4-s > 30
    [ask n-of 30 a-4-s
      [set next-gen-tag 1]]
    [ask a-4-s
      [set next-gen-tag 1]]
    ;; a-5-s
     ifelse count a-5-s > 30
    [ask n-of 30 a-5-s
      [set next-gen-tag 1]]
    [ask a-4-s
      [set next-gen-tag 1]]
    ;; a-1-r
    ifelse count a-1-r > 30
    [ask n-of 30 a-1-r
      [set next-gen-tag 1]]
    [ask a-1-r
      [set next-gen-tag 1]]
    ifelse count a-2-r > 30
    [ask n-of 30 a-2-r
      [set next-gen-tag 1]]
    [ask a-2-r
      [set next-gen-tag 1]]
    ifelse count a-3-r > 30
    [ask n-of 30 a-3-r
      [set next-gen-tag 1]]
    [ask a-3-s
      [set next-gen-tag 1]]
    ifelse count a-4-r > 30
    [ask n-of 30 a-4-r
      [set next-gen-tag 1]]
    [ask a-4-r
      [set next-gen-tag 1]]
     ifelse count a-5-r > 30
    [ask n-of 30 a-5-r
      [set next-gen-tag 1]]
    [ask a-5-r
      [set next-gen-tag 1]]

    ;;;;;;;;;;;;;;;
    ;; substrates
    ;;;;;;;;;;;;;;;
    ifelse count substrate > 300
    [ask n-of 30 substrate
      [set next-gen-tag 1]]
    [ask substrate
      [set next-gen-tag 1]]

    ;; mod-a-1-n
    ifelse count mod-a-1-n > 300
    [ask n-of 30 mod-a-1-n
      [set next-gen-tag 1]]
    [ask mod-a-1-n
      [set next-gen-tag 1]]


    ifelse count mod-a-2-n > 300
    [ask n-of 30 mod-a-2-n
      [set next-gen-tag 1]]
    [ask mod-a-2-n
      [set next-gen-tag 1]]


    ifelse count mod-a-3-n > 300
    [ask n-of 30 mod-a-3-n
      [set next-gen-tag 1]]
    [ask mod-a-3-n
      [set next-gen-tag 1]]


    ifelse count mod-a-4-n > 300
    [ask n-of 30 mod-a-4-n
      [set next-gen-tag 1]]
    [ask mod-a-4-n
      [set next-gen-tag 1]]

    ifelse count substrate-s > 300
    [ask n-of 30 substrate-s
      [set next-gen-tag 1]]
    [ask substrate-s
      [set next-gen-tag 1]]

    ifelse count mod-a-1-s > 300
    [ask n-of 30 mod-a-1-s
      [set next-gen-tag 1]]
    [ask mod-a-1-s
      [set next-gen-tag 1]]


    ifelse count mod-a-2-s > 300
    [ask n-of 30 mod-a-2-s
      [set next-gen-tag 1]]
    [ask mod-a-2-s
      [set next-gen-tag 1]]

    ifelse count mod-a-3-s > 300
    [ask n-of 30 mod-a-3-s
      [set next-gen-tag 1]]
    [ask mod-a-3-s
      [set next-gen-tag 1]]


    ifelse count mod-a-4-s > 300
    [ask n-of 30 mod-a-4-s
      [set next-gen-tag 1]]
    [ask mod-a-4-s
      [set next-gen-tag 1]]


    ask turtles
    [if not (next-gen-tag = 1)
      [die]]
    ]
  ;; Any changes within global gamma-lambda must happen here as well
  birth-mut-global
  birth-mut-counter-proc

  implicit-interaction-here-sub
  implicit-interaction-here-a-1-n
  implicit-interaction-here-a-2-n
  implicit-interaction-here-a-3-n
  implicit-interaction-here-a-4-n
  implicit-interaction-here-a-5-n
  implicit-interaction-here-mod-a-1-n
  implicit-interaction-here-mod-a-2-n
  implicit-interaction-here-mod-a-3-n
  implicit-interaction-here-mod-a-4-n
  implicit-interaction-here-a-1-r
  implicit-interaction-here-a-2-r
  implicit-interaction-here-a-3-r
  implicit-interaction-here-a-4-r
  implicit-interaction-here-a-5-r
  implicit-interaction-here-a-1-s
  implicit-interaction-here-a-2-s
  implicit-interaction-here-a-3-s
  implicit-interaction-here-a-4-s
  implicit-interaction-here-a-5-s
  implicit-interaction-here-sub-s


  implicit-interaction-gone-sub
  implicit-interaction-gone-a-1-n
  implicit-interaction-gone-a-2-n
  implicit-interaction-gone-a-3-n
  implicit-interaction-gone-a-4-n
  implicit-interaction-gone-a-5-n
  implicit-interaction-gone-a-1-s
  implicit-interaction-gone-a-2-s
  implicit-interaction-gone-a-3-s
  implicit-interaction-gone-a-4-s
  implicit-interaction-gone-a-5-s
  implicit-interaction-gone-a-1-r
  implicit-interaction-gone-a-2-r
  implicit-interaction-gone-a-3-r
  implicit-interaction-gone-a-4-r
  implicit-interaction-gone-a-5-r
  implicit-interaction-gone-mod-a-1-n
  implicit-interaction-gone-mod-a-2-n
  implicit-interaction-gone-mod-a-3-n
  implicit-interaction-gone-mod-a-4-n
  implicit-interaction-gone-sub-s



  tick-advance 1
  display



end

to move
  if patch-ahead (speed * tick-delta) != patch-here
    [ set last-collision nobody ]
  jump (speed * tick-delta)
  set heading random-float 360
end
@#$#@#$#@
GRAPHICS-WINDOW
216
10
628
443
100
100
2.0
1
10
1
1
1
0
1
1
1
-100
100
-100
100
0
0
1
ticks
100.0

SLIDER
37
235
209
268
temperature
temperature
0
2000
200
100
1
NIL
HORIZONTAL

BUTTON
31
47
97
80
NIL
setup
NIL
1
T
OBSERVER
NIL
Z
NIL
NIL
1

BUTTON
31
79
97
112
NIL
go
T
1
T
OBSERVER
NIL
X
NIL
NIL
1

MONITOR
37
311
130
356
NIL
performance
17
1
11

MONITOR
37
268
130
313
NIL
count mod-a-1-n
17
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

lily pad
false
0
Polygon -7500403 true true 148 151 137 37 119 25 111 36 78 54 40 99 30 137 32 175 56 223 87 251 137 275 157 275 213 250 239 221 257 178 262 137 244 91 210 53 172 37 160 22 154 36
Line -13840069 false 154 151 207 97
Circle -13840069 false false 133 148 26
Line -13840069 false 52 122 134 157
Line -13840069 false 133 171 89 196
Line -13840069 false 147 193 147 254
Line -13840069 false 157 171 205 233
Line -13840069 false 161 161 204 163
Line -13840069 false 141 149 111 72

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

starry
false
0
Circle -7500403 false true 45 45 210
Polygon -7500403 true true 96 225 150 60 206 224 63 120 236 120
Polygon -7500403 true true 120 120 195 120 180 180 180 185 113 183

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

ufo top
false
0
Circle -2064490 true false 15 15 270
Circle -7500403 true true 75 75 150
Circle -16777216 false false 75 75 150
Circle -7500403 true true 60 60 30
Circle -7500403 true true 135 30 30
Circle -7500403 true true 210 60 30
Circle -7500403 true true 240 135 30
Circle -7500403 true true 210 210 30
Circle -7500403 true true 135 240 30
Circle -7500403 true true 60 210 30
Circle -7500403 true true 30 135 30
Circle -16777216 false false 30 135 30
Circle -16777216 false false 60 210 30
Circle -16777216 false false 135 240 30
Circle -16777216 false false 210 210 30
Circle -16777216 false false 240 135 30
Circle -16777216 false false 210 60 30
Circle -16777216 false false 135 30 30
Circle -16777216 false false 60 60 30

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.3.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="gen" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1600"/>
    <metric>ticks</metric>
    <metric>performance</metric>
  </experiment>
  <experiment name="BARE" repetitions="20" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>ticks</metric>
    <metric>performance</metric>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
