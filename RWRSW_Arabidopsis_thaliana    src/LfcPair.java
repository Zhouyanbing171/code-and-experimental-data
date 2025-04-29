import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;

public class LfcPair {
    public static HashMap<String, HashSet<String>> getLfc(String gene_type, ArrayList<Annotation> annotationList, ArrayList<Term> obo_terms) throws IOException {
        HashMap<String, HashSet<String>> ecNumGeneHashMap;

        // 基因证据ec编号 哈希表 初始化
        if (gene_type.equals("Arabidopsis_thaliana")) {
            ecNumGeneHashMap = Reader.getArabidopsisEcNumGeneHashMap(Calculator.getFilePaths(gene_type));
        } else {
            ecNumGeneHashMap = new HashMap<>();
        }
        // 将 obo_terms 转换为 Map，便于快速查找
        Map<String, String> goTermNamespaceMap = new HashMap<>();
        for (Term term : obo_terms) {
            goTermNamespaceMap.put(term.id, term.namespace);
        }
        // 分别遍历ecHash表和annotationList，将annotationList中的标识符ID和基因与ec表的标识符ID和基因校对
        // 对同一基因不同标识符ID的情况，以annotationList中的标识符ID为准
        HashMap<String, HashSet<String>> ecNumNewGeneHashMap = new HashMap<>();
        for (String ecNum : ecNumGeneHashMap.keySet()) {
            // 确保每个 ecNum 都有一个 HashSet 实例
            ecNumNewGeneHashMap.putIfAbsent(ecNum, new HashSet<>());
            for (String geneNameId : ecNumGeneHashMap.get(ecNum)) {
                String[] geneNameIdArr = geneNameId.split("_");
                String gene = geneNameIdArr[0];
                HashSet<String> annotationGoList = new HashSet<>();
                boolean isFound = false;
                // 分别遍历 annotationList 和 obo_terms 表，通过geneNameIdArr查找annotationList里对应的gene，用geneNameIdArr[0]找到的annotationList里对应的geneName就指定为gene，
                // 否则用geneNameIdArr[1]去找annotationList里对应的geneID，再据此找到的geneName就指定为gene，
                // 再查找这个gene对应的所有的goID总共是否有3种及以上的namespace，若没有则跳过，
                // 若有则添加到ecNumNewGeneHashMap中，名字以annotationList中的geneName为准
                for (Annotation annotation : annotationList) {
                    String geneName = annotation.getGeneName();
                    String geneID = annotation.getGeneID();
                    String goID = annotation.getGoID();
                    if (geneNameIdArr[0].equals(geneName)) {
                        annotationGoList.add(goID);
                        continue;
                    }
                    if (geneNameIdArr[1].equals(geneID)) {
                        annotationGoList.add(goID);
                        if (!isFound) {
                            if (ecNumNewGeneHashMap.get(ecNum).contains(geneName)) {
                                gene = geneID;
                                isFound = true;
                                continue;
                            }
                            gene = geneName;
                            isFound = true;
                        }
                    }
                }
                // 确保annotationGoList不为空
                if (!annotationGoList.isEmpty()) {
                    // 找到了对应的gene，开始查找goID
                    HashSet<String> goNamespaces = new HashSet<>();
                    for (String go : annotationGoList) {
                        String goNamespace = goTermNamespaceMap.get(go);
                        if (goNamespace != null) {
                            goNamespaces.add(goNamespace);
                        }
                    }
                    if (goNamespaces.size() >= 3) {
                        ecNumNewGeneHashMap.get(ecNum).add(gene);
                    }
                }
            }
        }
        ecNumGeneHashMap = ecNumNewGeneHashMap;

        // 过滤掉只包含一个基因的EC编号
        ecNumGeneHashMap.values().removeIf(entry -> entry.size() <= 1);

        // 创建一个HashSet来存储所有唯一的基因
        HashSet<String> uniqueGenes = new HashSet<>();
        ecNumGeneHashMap.values().forEach(uniqueGenes::addAll);

        // 转换为List以便后续处理
        List<String> uniqueGenesList = new ArrayList<>(uniqueGenes);

        // 确保文件路径存在
        final Path path = Paths.get("temp/" + gene_type + "/lfcPair.txt");
        Calculator.createFileIfNotExists(path);
        // 使用 TRUNCATE_EXISTING 确保覆盖原有内容
        Files.write(path, new byte[0], StandardOpenOption.TRUNCATE_EXISTING); // 先清空文件内容

        try (BufferedWriter writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            StringBuilder stringBuilder = new StringBuilder();
            int pairCount = 0; // 计数器，用于控制写入频率

            // 使用普通循环创建基因对组合并写入文件
            for (int i = 0; i < uniqueGenesList.size(); i++) {
                for (int j = i; j < uniqueGenesList.size(); j++) {
                    // 生成基因对
                    stringBuilder.append(uniqueGenesList.get(i))
                            .append("!")
                            .append(uniqueGenesList.get(j))
                            .append(System.lineSeparator());
                    pairCount++;

                    // 每1000对写入一次
                    if (pairCount % 1000 == 0) {
                        writer.write(stringBuilder.toString());
                        stringBuilder.setLength(0); // 清空StringBuilder
                        writer.flush(); // 刷新缓冲区
                    }
                }
            }

            // 写入剩余的基因对
            if (stringBuilder.length() > 0) {
                writer.write(stringBuilder.toString());
                writer.flush(); // 刷新缓冲区
            }
        }
        return ecNumGeneHashMap;
    }
}
