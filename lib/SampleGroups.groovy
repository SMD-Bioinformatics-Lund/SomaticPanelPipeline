class SampleGroups {
    private static final log = org.slf4j.LoggerFactory.getLogger(SampleGroups)

    static boolean skip(Object group, String reason) {
        log.warn("Skipping sample group '${group}': ${reason}")
        false
    }

    static boolean isTumor(Object meta) {
        meta instanceof Map && meta.type in ['T', 'tumor']
    }

    static boolean isNormal(Object meta) {
        meta instanceof Map && meta.type in ['N', 'normal']
    }

    static List samples(Object meta) {
        def samples = meta instanceof Map ? [meta] : meta
        if (!(samples instanceof List) || !samples || !samples.every { it instanceof Map }) {
            return []
        }
        samples
    }

    static List validateMetadata(Object group, Object meta) {
        def members = samples(meta)
        if (!members) {
            skip(group, 'metadata must be a map or a nonempty list of maps')
            return []
        }
        def ids = members.collect { it.id }
        if (ids.any { !it } || ids.unique(false).size() != ids.size() ||
            members.any { it.group != group || (!isTumor(it) && !isNormal(it)) } ||
            members.size() > 2 ||
            (members.size() == 2 && (members.count { isTumor(it) } != 1 || members.count { isNormal(it) } != 1))) {
            skip(group, "expected one sample or one tumor/normal pair with distinct IDs and matching metadata groups; found ${ids}")
            return []
        }
        members
    }

    static List validate(List record) {
        if (!record || record.size() < 2) {
            skip(record ? record[0] : null, 'missing group or metadata columns')
            return null
        }
        def samples = validateMetadata(record[0], record[1])
        if (!samples) {
            return null
        }
        if (record.drop(2).any { !(it instanceof List) || it.size() != samples.size() }) {
            skip(record[0], 'unaligned sample/file lists')
            return null
        }
        record
    }

    static boolean hasTumor(Object meta) {
        samples(meta).any { isTumor(it) }
    }

    static boolean matchesMember(Object group, Object meta, Object caseMeta) {
        if (!(meta instanceof Map)) {
            return skip(group, 'per-sample metadata must be a map')
        }
        def members = validateMetadata(group, caseMeta)
        if (!members) {
            return false
        }
        def member = members.find { it.id == meta.id }
        if (meta.group != group || !member ||
            isTumor(meta) != isTumor(member) || isNormal(meta) != isNormal(member)) {
            return skip(group, "sample '${meta.id}' does not belong to its germline VCF metadata")
        }
        true
    }
}
