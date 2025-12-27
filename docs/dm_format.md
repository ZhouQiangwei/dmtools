# dm On-Disk Format

## 概述

`dm` 是由 dmtools 生成的二进制甲基化覆盖度文件，沿用 bigWig/binaMeth 的封装：按染色体分块写入压缩的数据块，并辅以 R-tree 索引以支持区间查询。

设计目标：

- 支持远程和本地随机访问（依赖块级索引）。
- 跨平台稳定：所有字段均为固定宽度并以 little-endian 编码，不依赖 `struct` 对齐或 `sizeof(struct)`。
- 通过 `magic`、`version/type` 掩码与索引自检快速发现损坏。

## 文件头（Header）

所有整数以 little-endian 存储，文件开头与结尾各写一份 `BIGWIG_MAGIC` (`0x888FFC26`) 作为哨兵。

`binaMethHdr_t.version` 在内存中固定为 `4`，但 **磁盘第 0x04 偏移处写入的是 `uint16_t version/type_mask = BM_MAGIC | BM_*`**。默认仅有 `BM_MAGIC` (`0x8000`) 置位；命令行 `-C/-S/--Cx/--Id/-E` 会额外置位 `BM_COVER`、`BM_STRAND`、`BM_CONTEXT`、`BM_ID`、`BM_END`，并影响记录布局与 `fieldCount` 计算。

| 偏移 | 类型       | 字段                 | 说明 |
| ---- | ---------- | -------------------- | ---- |
| 0x00 | `uint32_t` | `magic`              | 必为 `BIGWIG_MAGIC`。文件末尾也会写一份作为尾标。 |
| 0x04 | `uint16_t` | `version/type_mask`  | `BM_MAGIC` 与可选 `BM_*` 位的按位或结果，用于描述数据记录包含的字段。 |
| 0x06 | `uint16_t` | `nLevels`            | 预留的缩放层级数量（默认 0，或由 `--zl` 设置）。 |
| 0x08 | `uint64_t` | `ctOffset`           | 染色体名称/长度树的起始偏移。 |
| 0x10 | `uint64_t` | `dataOffset`         | 数据区起始偏移；该位置首先存放 `uint64_t nBlocks`。 |
| 0x18 | `uint64_t` | `indexOffset`        | 主索引（R-tree）起始偏移。 |
| 0x20 | `uint16_t` | `fieldCount`         | 数据记录包含的字段数量，依据 `BM_*` 掩码计算。 |
| 0x22 | `uint16_t` | `definedFieldCount`  | 与 `fieldCount` 相同。 |
| 0x24 | `uint64_t` | `sqlOffset`          | 未使用，默认为 0。 |
| 0x2c | `uint64_t` | `summaryOffset`      | 指向摘要块；摘要为 40 字节：`uint64_t nBasesCovered`、`double minVal`、`double maxVal`、`double sumData`、`double sumSquared`。 |
| 0x34 | `uint32_t` | `bufSize`            | 数据块压缩缓冲区大小（默认 32768）。非零表示数据块使用 zlib 压缩。 |
| 0x38 | `uint64_t` | `extensionOffset`    | 可选扩展区起始偏移（例如记录写入参数）。 |
| 0x40 | 24×`nLevels` 字节 | Zoom headers | 每级 24 字节：`uint32_t level`、`uint32_t padding`、`uint64_t dataOffset`、`uint64_t indexOffset`。当前工具通常为 0。 |

### 染色体列表（`ctOffset`）
使用 `CIRTREE_MAGIC` (`0x78CA8C91`) 的树结构保存 contig 名称和长度。键为染色体字符串（不以 `\0` 结尾）；值为 `uint32_t` 长度。树节点以固定宽度字段写入，不依赖结构体对齐。

## 基础类型约定

- 所有 on-disk 字段均使用固定宽度类型：`uint8_t`、`uint16_t`、`uint32_t`、`uint64_t`、`float`、`double`。
- 禁止直接 `fwrite` 整个内存结构体；写入顺序按字段逐个拷贝（参见 `dmWrite.c` 中的 `memcpy` 序列）。
- 编码固定为 little-endian；其他字节序需在读取时自行转换。

## 主体数据区（Data Records）

数据区位于 `dataOffset`，开头是 `uint64_t nBlocks`，随后为 `nBlocks` 个数据块。每个数据块：

1. **块头（24 字节）**
   - `uint32_t tid`：染色体在染色体列表中的索引。
   - `uint32_t start`：块内最小起始坐标（0-based，包含）。
   - `uint32_t end`：块内最大结束坐标（不包含）。
   - `uint32_t step`：步长（当前写入固定为 0，表示可变步长）。
   - `uint32_t span`：跨度（当前写入固定为 0）。
   - `uint8_t ltype`：块类型；bam2dm 始终写入 `1`（bedGraph 型可变步长记录）。
   - `uint8_t padding`：保留，固定 0。
   - `uint16_t nItems`：块内记录条数。

2. **记录列表**（`nItems` 条）。字段顺序由 `version/type_mask` 中的 `BM_*` 位决定：
   - `uint32_t start`（必选；**当前 dmtools 写入 1-based 坐标**）。
   - `uint32_t end`（当 `BM_END` 置位时写入；end 为半开区间的右端点）。
   - `float value`（必选；甲基化比例）。
   - `uint16_t coverage`（当 `BM_COVER` 置位时写入；测序覆盖度）。
   - `uint8_t strand`（当 `BM_STRAND` 置位时写入；`'+'`/`'-'` 以 `0/1` 表示）。
   - `uint8_t context`（当 `BM_CONTEXT` 置位时写入；碱基上下文编码）。
   - `uint32_t entryid`（当 `BM_ID` 置位时写入；单个 4 字节数值 ID）。

当 `BM_ID` 置位时，dmtools 还会生成一个同名侧写文件 `<dm>.idmap.tsv`，用于保存数值 ID 与原始细胞条形码/名称的映射。

记录在文件内按 `tid`、`start` 的写入顺序保存；索引假定这种单调性以支持区间查询。`bufSize` 非 0 时，块（含块头+记录）经 zlib 压缩写入，索引保存压缩后块的偏移与大小。

## 扩展区（`extensionOffset`）

若 `extensionOffset` 非 0，dmtools 会写入一个可选扩展结构记录写入参数，便于重现和排障。当前扩展版本包含：

- `uint32_t magic`：`0x44574d50`（"DWMP"）
- `uint16_t version`：`1`
- `uint16_t size`：结构体字节数
- `uint32_t bufSize`：写入时使用的压缩缓冲区大小
- `uint32_t blockSize`：索引节点最大子节点数

## 索引结构（Index）

`indexOffset` 指向数据块的 R-tree 索引，头部格式：

- `uint32_t magic`：`IDX_MAGIC` (`0x2468ACE0`)。
- `uint32_t blockSize`：索引节点最大子节点数（默认 256，可通过 CLI 覆盖）。
- `uint64_t nBlocks`：数据块数量，应与 `nBlocks` 计数一致。
- `uint32_t startTid`、`uint32_t startPos`、`uint32_t endTid`、`uint32_t endPos`：根节点覆盖的范围。
- `uint64_t indexSize`：索引总字节数（包含头部和所有节点）。
- `uint32_t itemsPerBlock`：固定写入 1（与 bigWig 一致）。
- `uint32_t padding`：保留 0。
- `uint64_t rootOffset`：根节点起始偏移（紧随头部，8 字节对齐）。

节点格式：

- `uint8_t isLeaf`、`uint8_t padding`、`uint16_t nChildren`。
- 对每个子节点：
  - `uint32_t chrIdxStart`、`uint32_t baseStart`、`uint32_t chrIdxEnd`、`uint32_t baseEnd`。
  - 若叶子：`uint64_t dataOffset`、`uint64_t dataSize`（指向压缩块）。
  - 若非叶：`uint64_t childOffset` 指向子节点。

索引节点按 `chrIdx/start` 升序组织，允许根据 `tid/start-end` 区间快速定位块，实现随机 seek。

## 完整性与校验

- 读取时应验证：
  - 头部和尾部的 `magic` 均为 `BIGWIG_MAGIC`。
  - `version/type_mask` 至少包含 `BM_MAGIC`；若缺失或与字段数量不一致应视为损坏。
  - `dataOffset`、`indexOffset`、`ctOffset`、`summaryOffset` 均落在文件范围内。
  - `nBlocks` 与索引头的 `nBlocks` 一致，索引叶子中的 `(offset, size)` 不越界且按文件顺序非递减。
- 索引假设数据块按 `tid`、`start` 递增；违背该假设可能导致迭代失败或空结果。

## 向后兼容策略

- `version/type_mask` 中的 `BM_*` 位描述记录布局；新增字段应分配新的 `BM_*` 位并同步调整 `fieldCount`/写入顺序。
- 只要 `BM_MAGIC` 保持、索引与数据偏移合法，旧版本 reader 可以忽略未知 `BM_*` 位或拒绝解析，以避免误读数据偏移。
- 当前自检模式（`bam2dm --check` / `--validate-output`）视 `magic` 不匹配、偏移越界或索引非单调为损坏。
