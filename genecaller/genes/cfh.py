"""CFH gene with enhanced interpretation for CFHR deletions."""

from typing import Dict, Any
from genecaller.gene import Gene


class CFH(Gene):
    """CFH gene with enhanced interpretation for CFHR deletions."""

    def prepare_output(self) -> Dict[str, Any]:
        result_data = super().prepare_output()

        interpretation_lines = []
        copy_numbers = result_data["Copy numbers"]

        # Track CNV events
        cnv_events = []
        for r in self.cn_regions:
            region_name, coords = r
            cn = copy_numbers.get(region_name, 2)

            if cn == 0:
                event_type = "Homozygous deletion"
                cnv_events.append((region_name, event_type, coords))
                interpretation_lines.append(
                    f"{region_name}: {event_type} (CN=0, {coords})"
                )
            elif cn == 1:
                event_type = "Heterozygous deletion"
                cnv_events.append((region_name, event_type, coords))
                interpretation_lines.append(
                    f"{region_name}: {event_type} (CN=1, {coords})"
                )
            elif cn == 2:
                interpretation_lines.append(f"{region_name}: Normal copy number (CN=2)")
            elif cn == 3:
                event_type = "Heterozygous duplication"
                cnv_events.append((region_name, event_type, coords))
                interpretation_lines.append(
                    f"{region_name}: {event_type} (CN=3, {coords})"
                )
            else:
                event_type = f"Copy number variation"
                cnv_events.append((region_name, event_type, coords))
                interpretation_lines.append(
                    f"{region_name}: {event_type} (CN={cn}, {coords})"
                )

        # Add summary if CNV events detected
        if cnv_events:
            interpretation_lines.append(f"CNV events detected: {len(cnv_events)}")

        result_data["cn_interpretation"] = "\n".join(interpretation_lines)

        return result_data
