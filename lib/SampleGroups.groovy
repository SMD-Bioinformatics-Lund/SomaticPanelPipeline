class SampleGroups {
    static boolean isTumor(Map meta) {
        meta.type in ['T', 'tumor']
    }

    static boolean isNormal(Map meta) {
        meta.type in ['N', 'normal']
    }

    static List samples(Object meta) {
        def samples = meta instanceof Map ? [meta] : meta
        if (!(samples instanceof List) || !samples || !samples.every { it instanceof Map }) {
            throw new IllegalArgumentException('Sample metadata must be a map or a nonempty list of maps')
        }
        samples
    }

    static List validateMetadata(Object group, Object meta) {
        def members = samples(meta)
        def ids = members.collect { it.id }
        if (ids.any { !it } || ids.unique(false).size() != ids.size() ||
            members.any { it.group != group || (!isTumor(it) && !isNormal(it)) } ||
            members.size() > 2 ||
            (members.size() == 2 && (members.count { isTumor(it) } != 1 || members.count { isNormal(it) } != 1))) {
            throw new IllegalArgumentException("Invalid sample group '${group}': expected one sample or one tumor/normal pair with distinct IDs and matching metadata groups; found ${ids}")
        }
        members
    }

    static List validate(List record) {
        def samples = validateMetadata(record[0], record[1])
        if (record.drop(2).any { !(it instanceof List) || it.size() != samples.size() }) {
            throw new IllegalArgumentException("Unaligned sample/file lists in group '${record[0]}'")
        }
        record
    }

    static boolean hasTumor(Object meta) {
        samples(meta).any { isTumor(it) }
    }

    static void requireMember(Object group, Map meta, Object caseMeta) {
        def members = validateMetadata(group, caseMeta)
        def member = members.find { it.id == meta.id }
        if (meta.group != group || !member ||
            isTumor(meta) != isTumor(member) || isNormal(meta) != isNormal(member)) {
            throw new IllegalArgumentException("Sample '${meta.id}' does not belong to germline VCF group '${group}'")
        }
    }
}
