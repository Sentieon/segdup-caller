"""CFH gene with enhanced interpretation for CFHR deletions."""

from typing import Dict, Any
from genecaller.gene import Gene


class CFH(Gene):
    """CFH gene with enhanced interpretation for CFHR deletions."""

    def prepare_output(self) -> Dict[str, Any]:
        result_data = super().prepare_output()

        interpretation_lines = []
        copy_numbers = result_data["Copy numbers"]

        # Build region_info from all_vars["cns"] which contains cn_regions from config
        region_info = {}
        for region_data in self.all_vars["cns"].values():
            region_name = region_data["name"]
            coords = self.cn_regions[region_name]

            # Parse coordinates to calculate size
            _, positions = coords.split(":")
            start, end = map(int, positions.split("-"))
            size_bp = end - start

            region_info[region_name] = (coords, size_bp)

        # Track CNV events
        cnv_events = []

        for region_name, (coords, size_bp) in region_info.items():
            cn = copy_numbers.get(region_name, 2)
            size_kb = size_bp / 1000

            if cn == 0:
                event_type = "Homozygous deletion"
                cnv_events.append((region_name, event_type, coords, size_kb))
                interpretation_lines.append(
                    f"{region_name}: {event_type} (CN=0, {coords}, {size_kb:.1f} kb)"
                )
            elif cn == 1:
                event_type = "Heterozygous deletion"
                cnv_events.append((region_name, event_type, coords, size_kb))
                interpretation_lines.append(
                    f"{region_name}: {event_type} (CN=1, {coords}, {size_kb:.1f} kb)"
                )
            elif cn == 2:
                interpretation_lines.append(f"{region_name}: Normal copy number (CN=2)")
            elif cn == 3:
                event_type = "Heterozygous duplication"
                cnv_events.append((region_name, event_type, coords, size_kb))
                interpretation_lines.append(
                    f"{region_name}: {event_type} (CN=3, {coords}, {size_kb:.1f} kb)"
                )
            else:
                event_type = f"Copy number variation"
                cnv_events.append((region_name, event_type, coords, size_kb))
                interpretation_lines.append(
                    f"{region_name}: {event_type} (CN={cn}, {coords}, {size_kb:.1f} kb)"
                )

        # Add summary if CNV events detected
        if cnv_events:
            interpretation_lines.append("")
            interpretation_lines.append(f"CNV events detected: {len(cnv_events)}")
            total_size_kb = sum(size for _, _, _, size in cnv_events)
            interpretation_lines.append(
                f"Total affected region: {total_size_kb:.1f} kb"
            )

        result_data["cn_interpretation"] = "\n".join(interpretation_lines)

        return result_data
