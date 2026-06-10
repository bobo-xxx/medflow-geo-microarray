# Core Connection

## Purpose

Every node lives in its own git repository and connects to the IRE core framework via a declarative `registry.yaml` in the core repo. This spec defines the linkage protocol — what the core needs from a node repo, and what the node author must provide.

## Requirements

### Requirement: Repository Structure

The node repo SHALL contain the package files in a `node/` subdirectory: `node/SKILL.md`, `node/env.yaml`, `node/scripts/`. This separates the node package (what the framework consumes) from project-level files (CLAUDE.md, docs/, openspec/, tests/, .github/). The version is encoded in the registry.yaml in the core, NOT in the repo directory name.

```
node-repo/
├── node/
│   ├── SKILL.md              # Agent contract (frontmatter + body)
│   ├── env.yaml              # Declarative conda environment
│   ├── scripts/              # Entry point + internal modules
│   └── references/           # Optional static assets
├── tests/                    # Test suite
├── docs/                     # Design docs, plans
├── openspec/                 # Project governance
├── .claude/                  # Claude Code configuration
├── .github/                  # CI workflows
├── CLAUDE.md                 # Developer guidance
└── .gitignore
```

#### Scenario: Core clones a node

- **WHEN** the core's `ire sync` clones this repo
- **THEN** it SHALL checkout a pinned commit
- **AND** extract the `node/` subdirectory (if present) as the package root; if `node/` is absent, use the repo root
- **AND** place the package at `nodes/<name>@<version>/` relative to the core root
- **AND** verify sha256 integrity of the package files

### Requirement: Registry Entry

The node SHALL be declared in the core's `registry.yaml`:

```yaml
nodes:
  geo-microarray-processing:
    source: git
    url: https://github.com/ire-product/geo-microarray-processing.git
    versions:
      - version: "1.0.0"
        commit: "<full-git-sha>"
        sha256: "<sha256-of-package-files>"
```

#### Scenario: New version published

- **WHEN** a new version of this node is released
- **THEN** the node author SHALL tag the commit and provide the commit hash and sha256 to the core maintainer
- **AND** the core maintainer SHALL add a new entry under `versions:`

#### Scenario: Node author tests locally

- **WHEN** developing this node
- **THEN** the author SHALL symlink or copy the `node/` subdirectory to `nodes/<name>@<version>/` in a local core checkout
- **AND** the core's `NodeRegistry.discover()` SHALL index it for testing without `ire sync`

### Requirement: No Core Code in Node Repo

The node repo SHALL contain only node package files. Framework code, protocol documents, and other nodes live in their own repos.

#### Scenario: Node repo boundaries

- **WHEN** inspecting a node repo
- **THEN** it SHALL contain: `node/SKILL.md`, `node/env.yaml`, `node/scripts/`, `node/references/` (optional), `openspec/` (project governance)
- **AND** it SHALL NOT contain: `registry.yaml` (core), `specs/m2m/` (core), framework engine code, or other nodes

### Requirement: GitHub Token

The node author SHALL store the GitHub token in `.claude/settings.local.json` (never committed). It is available as `$GITHUB_TOKEN`.
