import org.yaml.snakeyaml.Yaml

def resolveFastqPath(String x) {
    if (x.startsWith('FASTQ') || x.startsWith('/')) {
        return x
    }
    return "FASTQ/${x}"
}

def loadAlignmentInputs(String fastqConfigPath, String sampleListPath, String referenceKey) {
    def yaml = new Yaml()
    def files = yaml.load(new File(fastqConfigPath).text) as Map

    files.keySet().toList().each { sample ->
        def val = files[sample] as Map
        def flist = (val['files'] as List).sort()
        def r1ByBase = [:]
        def r2ByBase = [:]
        flist.each { f ->
            def base = new File(f as String).name
            if (base.contains('_R1')) {
                def baseStem = base.replace('_R1.fastq.gz', '').replace('_R1.fq.gz', '')
                r1ByBase[baseStem] = f
            } else if (base.contains('_R2')) {
                def baseStem = base.replace('_R2.fastq.gz', '').replace('_R2.fq.gz', '')
                r2ByBase[baseStem] = f
            }
        }
        def paired = []
        (r1ByBase.keySet() + r2ByBase.keySet()).unique().sort().each { baseStem ->
            def r1 = r1ByBase[baseStem]
            def r2 = r2ByBase[baseStem]
            if (r1 == null) {
                System.err.println("WARNING: unpaired R2 (no R1): ${r2}")
                return
            }
            if (r2 == null) {
                System.err.println("WARNING: unpaired R1 (no R2): ${r1}")
                return
            }
            paired << [baseStem, r1, r2]
        }
        if (paired.isEmpty()) {
            System.err.println("WARNING: no paired FASTQs for sample ${sample}, skipping")
            files.remove(sample)
            return
        }
        def newVal = [:]
        paired.each { baseStem, r1, r2 ->
            def parts = (baseStem as String).split('_')
            def run, lane, index
            if (parts.size() >= 3) {
                run = parts[-3]
                lane = parts[-2]
                index = parts[-1]
            } else {
                run = 'run'
                lane = '0'
                index = '0'
            }
            def key = "${run}-${lane}-${index}"
            newVal[key] = [PU: "${run}-${lane}-${index}", files: [r1, r2]]
        }
        files[sample] = newVal
    }

    def samples = new File(sampleListPath).readLines().collect { it.trim() }.findAll { it }
    samples.each { sample ->
        if (!files.containsKey(sample)) {
            System.err.println("WARNING: sample ${sample} in sample.list not in fastq config, skipping")
        }
    }
    samples = samples.findAll { files.containsKey(it) }

    def lanes = []
    samples.each { sample ->
        def sampleMap = files[sample] as Map
        sampleMap.keySet().sort().each { key ->
            def entry = sampleMap[key] as Map
            def puParts = (key as String).split('-', 3)
            def run = puParts[0]
            def lane = puParts.length > 1 ? puParts[1] : '0'
            def index = puParts.length > 2 ? puParts[2] : '0'
            def r1 = resolveFastqPath(entry.files[0] as String)
            def r2 = resolveFastqPath(entry.files[1] as String)
            lanes << [sample, referenceKey, run, lane, index, r1, r2]
        }
    }

    return [samples: samples, lanes: lanes]
}
