from __future__ import annotations

from dataclasses import dataclass
from typing import List, Iterable, Dict, Tuple
from collections import defaultdict


# ======================
# 数据结构
# ======================

@dataclass(frozen=True)
class CrosslinkSite:
    chrom: str
    start: int
    end: int
    score: float
    strand: str


@dataclass
class ClusterRegion:
    chrom: str
    start: int
    end: int
    strand: str
    sites: List[CrosslinkSite]


@dataclass
class BindingSite:
    chrom: str
    start: int
    end: int
    strand: str
    center: int
    support: int
    score: float


# ======================
# 参数配置（原则驱动）
# ======================

@dataclass(frozen=True)
class BindingConfig:
    extend: int = 4
    min_cluster_width: int = 3
    min_support: int = 3

    @property
    def target_width(self) -> int:
        return 2 * self.extend + 1

    @property
    def max_gap(self) -> int:
        # 原则：resize 后不重叠
        return self.target_width - 1


# ======================
# 解析
# ======================

def parse_pureclip(path: str) -> List[CrosslinkSite]:
    sites: List[CrosslinkSite] = []

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")

            sites.append(
                CrosslinkSite(
                    chrom=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    score=float(fields[4]),
                    strand=fields[5],
                )
            )

    return sites


# ======================
# Step 1: 聚类
# ======================

def cluster_sites(
    sites: Iterable[CrosslinkSite],
    config: BindingConfig
) -> List[ClusterRegion]:

    grouped: Dict[Tuple[str, str], List[CrosslinkSite]] = defaultdict(list)

    for s in sites:
        grouped[(s.chrom, s.strand)].append(s)

    clusters: List[ClusterRegion] = []

    for (chrom, strand), group in grouped.items():
        group.sort(key=lambda x: x.start)

        current_sites = [group[0]]
        start = group[0].start
        end = group[0].end

        for site in group[1:]:
            gap = site.start - end

            if gap <= config.max_gap:
                current_sites.append(site)
                end = max(end, site.end)
            else:
                clusters.append(
                    ClusterRegion(chrom, start, end, strand, current_sites)
                )
                current_sites = [site]
                start = site.start
                end = site.end

        clusters.append(
            ClusterRegion(chrom, start, end, strand, current_sites)
        )

    return clusters


# ======================
# Step 2: 过滤短 cluster
# ======================

def filter_clusters(
    clusters: List[ClusterRegion],
    config: BindingConfig
) -> List[ClusterRegion]:

    return [
        c for c in clusters
        if (c.end - c.start) >= config.min_cluster_width
    ]


# ======================
# 工具函数
# ======================

def center_window(center: int, width: int) -> Tuple[int, int]:
    half = width // 2
    return center - half, center + half + 1


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return not (a_end <= b_start or a_start >= b_end)


def overlaps_any(
    start: int,
    end: int,
    intervals: List[Tuple[int, int]]
) -> bool:
    for s, e in intervals:
        if overlaps(start, end, s, e):
            return True
    return False


def count_support(
    sites: List[CrosslinkSite],
    start: int,
    end: int
) -> int:
    return sum(1 for s in sites if start <= s.start < end)


# ======================
# Step 3: 统一宽度
# ======================

def generate_binding_sites(
    cluster: ClusterRegion,
    config: BindingConfig
) -> List[BindingSite]:

    width = cluster.end - cluster.start
    target = config.target_width

    # 按信号排序（高 → 低）
    sites_sorted = sorted(cluster.sites, key=lambda x: x.score, reverse=True)

    results: List[BindingSite] = []
    used_intervals: List[Tuple[int, int]] = []

    # ---- long region → 拆分 ----
    if width > target:
        for s in sites_sorted:
            center = s.start
            start, end = center_window(center, target)

            # 限制 center 必须在 cluster 内
            if not (cluster.start <= center < cluster.end):
                continue

            if overlaps_any(start, end, used_intervals):
                continue

            support = count_support(cluster.sites, start, end)

            results.append(
                BindingSite(
                    chrom=cluster.chrom,
                    start=start,
                    end=end,
                    strand=cluster.strand,
                    center=center,
                    support=support,
                    score=s.score,
                )
            )

            used_intervals.append((start, end))

    # ---- short / equal → 单窗口 ----
    else:
        peak = max(cluster.sites, key=lambda x: x.score)
        center = peak.start
        start, end = center_window(center, target)

        support = count_support(cluster.sites, start, end)

        results.append(
            BindingSite(
                chrom=cluster.chrom,
                start=start,
                end=end,
                strand=cluster.strand,
                center=center,
                support=support,
                score=peak.score,
            )
        )

    return results


# ======================
# Step 4: 支持度过滤
# ======================

def filter_binding_sites(
    sites: List[BindingSite],
    config: BindingConfig
) -> List[BindingSite]:

    return [
        s for s in sites
        if s.support >= config.min_support
    ]


# ======================
# Pipeline
# ======================

def build_binding_sites(
    sites: List[CrosslinkSite],
    config: BindingConfig
) -> List[BindingSite]:

    clusters = cluster_sites(sites, config)
    clusters = filter_clusters(clusters, config)

    results: List[BindingSite] = []

    for c in clusters:
        results.extend(generate_binding_sites(c, config))

    return filter_binding_sites(results, config)


# ======================
# 输出（BED-like）
# ======================

def write_binding_sites(path: str, sites: List[BindingSite]) -> None:
    with open(path, "w") as f:
        for s in sites:
            f.write(
                f"{s.chrom}\t{s.start}\t{s.end}\t{s.score:.6f}\t{s.support}\t{s.strand}\n"
            )


# ======================
# CLI 示例入口（可删）
# ======================

def main(input_path: str, output_path: str):
    config = BindingConfig(
        extend=4,
        min_cluster_width=3,
        min_support=3
    )

    sites = parse_pureclip(input_path)
    binding_sites = build_binding_sites(sites, config)
    write_binding_sites(output_path, binding_sites)


if __name__ == "__main__":
    import sys
    main(sys.argv[1], sys.argv[2])