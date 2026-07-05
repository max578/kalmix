# API stability policy — kalmix

kalmix is pre-1.0. The exported surface is settling and follows a
semver-flavoured policy:

- **Patch releases** (`0.x.y` → `0.x.y+1`) never change the exported API.
- **Minor releases** (`0.x` → `0.x+1`) may introduce breaking changes. Every
  breaking change is announced first in `NEWS.md` under a *Breaking changes*
  heading, with the before and after spelled out (the 0.1.0 and 0.3.0 entries
  model the form), and lands only at a minor bump, never in a patch.
- **Deprecations before removals** where the old surface can coexist with the
  new one; a removal without a deprecation cycle happens only when the old
  behaviour was incorrect (a bug masquerading as an API).
- Everything documented with `@noRd`/`@keywords internal`, and every object
  accessed via `:::`, is internal and carries no stability promise.

From 1.0.0 the policy hardens to full semantic versioning: breaking changes
only at major bumps.
