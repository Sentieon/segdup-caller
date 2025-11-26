"""HBA gene with clinical interpretation for alpha-thalassemia."""

from typing import Dict, Any
from genecaller.gene import Gene


class HBA(Gene):
    """HBA gene with alpha-thalassemia clinical interpretation."""

    def prepare_output(self) -> Dict[str, Any]:
        result_data = super().prepare_output()

        interpretation_lines = []
        copy_numbers = result_data["Copy numbers"]

        region_info = {r[0]: r[1] for r in self.cn_regions}
        # Get copy numbers (use standard naming with Greek alpha)
        region_display_names = {
            "a3.7": "-α3.7",
            "a4.2": "-α4.2",
            "non-duplication": "non-duplication",
        }

        cn_a37 = copy_numbers.get("a3.7", 2)
        cn_a42 = copy_numbers.get("a4.2", 2)
        cn_nondup = copy_numbers.get("non-duplication", 2)

        # Track CNV events
        cnv_events = []

        # Report individual region copy numbers
        for region_name, cn in copy_numbers.items():
            display_name = region_display_names.get(region_name, region_name)
            coords = region_info.get(region_name, "Non-duplication region")
            if cn == 0:
                event_type = "Homozygous deletion"
                cnv_events.append((region_name, event_type, coords))
                interpretation_lines.append(
                    f"{display_name}: {event_type} (CN=0, {coords})"
                )
            elif cn == 1:
                event_type = "Heterozygous deletion"
                cnv_events.append((region_name, event_type, coords))
                interpretation_lines.append(
                    f"{display_name}: {event_type} (CN=1, {coords})"
                )
            elif cn == 2:
                interpretation_lines.append(
                    f"{display_name}: Normal copy number (CN=2)"
                )
            elif cn == 3:
                event_type = "Triplication"
                cnv_events.append((region_name, event_type, coords))
                interpretation_lines.append(
                    f"{display_name}: {event_type} (CN=3, {coords})"
                )
            else:
                event_type = f"Copy number variation"
                cnv_events.append((region_name, event_type, coords))
                interpretation_lines.append(
                    f"{display_name}: {event_type} (CN={cn}, {coords})"
                )

        # Add clinical interpretation based on CN pattern
        clinical_interpretation = self._get_clinical_interpretation(
            cn_a37, cn_a42, cn_nondup
        )
        interpretation_lines.append(
            f"Clinical interpretation: {clinical_interpretation}"
        )

        result_data["cn_interpretation"] = "\n".join(interpretation_lines)

        return result_data

    def _get_clinical_interpretation(
        self, cn_a37: int, cn_a42: int, cn_nondup: int
    ) -> str:
        """
        Determine clinical interpretation based on copy number pattern.

        Args:
            cn_a37: Copy number of -α3.7 deletion region
            cn_a42: Copy number of -α4.2 deletion region
            cn_nondup: Copy number of non-duplication region

        Returns:
            Clinical interpretation string with genotype and phenotype
        """
        # Check for large deletions first (based on non-duplication region)
        if cn_nondup == 0:
            if cn_a37 == 0 and cn_a42 == 0:
                return "Homozygous large deletion (--/--) - Hemoglobin Bart's hydrops fetalis (lethal)"
            else:
                return "Abnormal non-duplication CN=0 - Requires manual review"

        if cn_nondup == 1:
            if cn_a37 == 0 and cn_a42 == 0:
                return "Heterozygous large deletion, likely --SEA, --MED, --FIL, or --THAI (αα/--) - Silent carrier"
            elif cn_a37 == 1 and cn_a42 == 2:
                return "Compound -α3.7 with large deletion (-α3.7/--) - Hemoglobin H disease"
            elif cn_a37 == 2 and cn_a42 == 1:
                return "Compound -α4.2 with large deletion (-α4.2/--) - Hemoglobin H disease"
            else:
                return f"Abnormal non-duplication CN=1 with unusual pattern - Requires manual review"

        # Standard patterns (non-duplication CN=2)
        if cn_a37 == 2 and cn_a42 == 2:
            return "Normal/Reference (αα/αα) - Normal"

        if cn_a37 == 1 and cn_a42 == 2:
            return "Heterozygous -α3.7 deletion (αα/-α3.7) - Silent carrier"

        if cn_a37 == 2 and cn_a42 == 1:
            return "Heterozygous -α4.2 deletion (αα/-α4.2) - Silent carrier"

        if cn_a37 == 0 and cn_a42 == 2:
            return "Homozygous -α3.7 deletion (-α3.7/-α3.7) - Alpha-thalassemia trait"

        if cn_a37 == 2 and cn_a42 == 0:
            return "Homozygous -α4.2 deletion (-α4.2/-α4.2) - Alpha-thalassemia trait"

        if cn_a37 == 1 and cn_a42 == 1:
            return "Compound heterozygous (-α3.7/-α4.2) - Alpha-thalassemia trait"

        if cn_a37 == 3 and cn_a42 == 2:
            return "α-globin triplication involving -α3.7 region - Usually normal"

        if cn_a37 == 2 and cn_a42 == 3:
            return "α-globin triplication involving -α4.2 region - Usually normal"

        # Non-duplication CN > 2
        if cn_nondup > 2:
            return f"Complex rearrangement (non-dup CN={cn_nondup}) - Requires manual review"

        # Catch-all for other patterns
        return f"Unusual copy number pattern (-α3.7={cn_a37}, -α4.2={cn_a42}, non-dup={cn_nondup}) - Requires manual review"
