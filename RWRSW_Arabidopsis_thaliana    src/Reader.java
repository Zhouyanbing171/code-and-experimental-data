import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Reader {
    // 读取基因本体文件，获得基因本体数据数组
    public static ArrayList<Term> readOboFile(String obo) throws IOException {
        ArrayList<Term> termList = new ArrayList<>();
        String content = new String(Files.readAllBytes(Paths.get(obo)));
        Pattern termPattern = Pattern.compile("\\[Term][\\s\\S.+?]+?\\n\\n");
        Pattern idPattern = Pattern.compile("(?<=id: )GO:\\d{7}");
        Pattern namePattern = Pattern.compile("(?<=name: ).+?\\\\n");
        Pattern namespacePattern = Pattern.compile("(?<=namespace: ).+?\\n");
        Pattern isaPattern = Pattern.compile("(?<=is_a: )GO:\\d{7}");
        Pattern partofPattern = Pattern.compile("(?<=relationship: part_of )GO:\\d{7}");
        Pattern obsoletePattern = Pattern.compile("is_obsolete: true");
        Matcher m = termPattern.matcher(content);
        while (m.find()) {
            Term termNode = new Term();
            String term = m.group();

            Matcher idm = idPattern.matcher(term);
            if (idm.find()) {
                termNode.id = idm.group().replace(":", "").trim();
            }
            Matcher names = namePattern.matcher(term);
            if (names.find()) {
                termNode.name = names.group().replace(":", "").trim();
            }
            Matcher nsm = namespacePattern.matcher(term);
            if (nsm.find()) {
                termNode.namespace = nsm.group().replace(":", "").trim();
            }
            Matcher isam = isaPattern.matcher(term);
            while (isam.find()) {
                termNode.isID.add(isam.group().replace(":", "").trim());
            }
            Matcher part_ofm = partofPattern.matcher(term);
            while (part_ofm.find()) {
                termNode.partID.add(part_ofm.group().replace(":", "").trim());
            }
            termNode.isObsolete = obsoletePattern.matcher(term).find();
            termList.add(termNode);
        }

        return termList;
    }

    // 读取注释信息文件，获得注释信息数据数组
    public static ArrayList<Annotation> readGafFile(String gaf) throws IOException {
        ArrayList<Annotation> annotationList = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(gaf))) {
            String line;
            while ((line = br.readLine()) != null) {
                // 跳过以!开头的注释行
                if (line.startsWith("!")) {
                    continue;
                }
                String[] row = line.split("\t");
                if (row.length >= 12) {
                    Annotation annotationNode = new Annotation();
                    // 第1列：基因名字
                    annotationNode.geneName = row[1];
                    // 第2列：基因ID
                    annotationNode.geneID = row[2];
                    // 第3列：基因本体术语ID
                    annotationNode.goID = row[4].replace(":", "");
                    // 第4列：证据编码
                    annotationNode.evidenceCode = row[6];
                    // 第5列：基因本体术语所属分支
                    annotationNode.nameSpace = row[8];
                    // 第6列：基因本体术语的同义词
                    String synonymData = row[9];
                    annotationNode.synonym = new ArrayList<>();
                    if (synonymData != null && !synonymData.isEmpty()) { // 检查是否不是空
                        String[] synonyms = synonymData.split("\\|");
                        for (String synonym : synonyms) {
                            // 过滤重复同义词
                            if (synonym != null && !synonym.trim().isEmpty()) {
                                annotationNode.synonym.add(synonym);
                            }
                        }
                    }
                    annotationList.add(annotationNode);
                }
            }
        }
        return annotationList;
    }

    // 读取基因功能网络文件，获得基因网络节点数据数组
    public static ArrayList<FunctionNet> readNetFile(String net) throws IOException {
        ArrayList<FunctionNet> netList = new ArrayList<>();
        for (String line : Files.readAllLines(Paths.get(net))) {
            FunctionNet netNode = new FunctionNet();
            netNode.gene1 = line.split("\t")[0];
            netNode.gene2 = line.split("\t")[1];
            netNode.value = Double.parseDouble(line.split("\t")[2]);
            netList.add(netNode);
        }
        return netList;
    }

    // 读取基因相似度文件，获得基因相似度数据数组 哈希表
    public static HashMap<String, Double> getSimilarityResult(String gene_type, double alpha) throws IOException {
        HashMap<String, Double> similarityResult = new HashMap<>();
        for (String line : Files.readAllLines(Paths.get("result/" + gene_type + "/similarityResult" + alpha + ".txt"))) {
            String gene1 = line.split("\t")[0];
            String gene2 = line.split("\t")[1];
            Double value = Double.parseDouble(line.split("\t")[2]);
            similarityResult.put(gene1 + ":" + gene2, value);
        }
        return similarityResult;
    }

    // 读取Arabidopsis thaliana ec文件，获得ec数据
    public static HashMap<String, HashSet<String>> getArabidopsisEcNumGeneHashMap(ArrayList<String> filePaths) throws IOException {
        // 基因证据ec编号 哈希表 初始化
        HashMap<String, HashSet<String>> ecNumGeneHashMap = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePaths.get(3)))) {
            String line;
            // 如果有则跳过表头
            br.readLine();
            while ((line = br.readLine()) != null) {
                String[] columns = line.split("\t"); // 使用制表符作为分隔符读取文件
                if (columns.length < 8) {
                    continue;
                } // 确保行中有足够的列
                String ecNumber = columns[3].replace("EC-", "").trim(); // EC编号
                String geneName = columns[6].trim(); // 基因名字
                String geneId = columns[7].trim(); // 基因ID
                if (geneName.equals("unknown") || ecNumber.equals("-")) {
                    continue;
                }
                // 检查EC编号和基因名字是否为空
                if (!ecNumber.isEmpty() && !geneName.isEmpty()) {
                    ecNumGeneHashMap.putIfAbsent(ecNumber, new HashSet<>()); // 如果没有则创建新的集合
                    ecNumGeneHashMap.get(ecNumber).add(geneName + "_" + geneId); // 添加基因名字到相应的EC编号
                }
            }
        } catch (IOException e) {
            throw new IOException(e);
        }
        // 过滤掉只包含一个基因的EC编号
        ecNumGeneHashMap.values().removeIf(entry -> entry.size() <= 1);
        // 移除空字符串的键（如果存在）
        ecNumGeneHashMap.remove("");
        return ecNumGeneHashMap;
    }

}
