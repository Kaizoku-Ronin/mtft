# Execution notes

The first primary training run completed its training calculations but stopped
while serializing a NumPy Boolean in the ambiguity flag. No selection JSON was
written and no held-out evaluation ran. The flag and its gap were converted to
native Python scalar types; training was then rerun under the unchanged frozen
protocol. This was an output-serialization correction, not a model, grid,
threshold, or measurement-rule change.

The protocol's vertical-stacking notation was clarified by an appended note
before execution. The final protocol hash is recorded in all run artifacts.
